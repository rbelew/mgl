''' crawl_ADV: utilities to crawl ADVina FAAH files
	was part of get_ADInfo

@version 1.0
@date on 30 Aug 14
@author: rbelew@ucsd.edu
'''

import sys
# import getopt
import re
import glob
import os
import csv

import tempfile
import tarfile
import shutil
import json
import argparse
import datetime
from string import strip

import numpy as np

import pybel
ob = pybel.ob

from CADD.Raccoon2.HelperFunctionsN3P import pathToList
import config
import FAAHA

def getLines(filename, doStrip = False, removeEmpty=False):
	""" """
	f = open(filename, 'r')
	lines = f.readlines()
	f.close()
	if doStrip:
		lines = map(strip,lines)
	if removeEmpty:
		#lines = removeEmptyLines(lines)
		lines = [ l for l in lines if l.strip() ]
	return lines

# defaults
mode = "1" # ONLY pose for Vina
ebest = -999. # energy
eworst = -3.
cworst = 1. # cluster poses
cbest = 999.
pworst = 1. # cluster %
pbest = 100.
lworst = 0. # ligand efficiency
lbest = -99
DEBUG = False
pattern = "_VS.pdbqt" # default pattern for searching result files
#pattern = ".VS.pdbqt" # default pattern for searching result files
do_filter = False
recursive = False

## Shared with AutoDock

# AADict = {'ASP':'D', 'GLU':'E', 'LYS':'K', 'HIS':'H', 'ARG':'R',
#				'GLN':'Q', 'ASN':'N', 'SER':'S', 'ASX':'B', 'GLX':'Z',
#				'PHE':'F', 'TRP':'W', 'TYR':'Y',
#				'GLY':'G', 'ALA':'A', 'ILE':'I', 'LEU':'L', 'CYS':'C',
#				'MET':'M', 'THR':'T', 'VAL':'V', 'PRO':'P' }

# updated Dec'14 with additional codes
AADict = {'ASP':'D', 'GLU':'E', 'LYS':'K', 'HIS':'H', 'ARG':'R',
			   'GLN':'Q', 'ASN':'N', 'SER':'S', 'ASX':'B', 'GLX':'Z',
			   'PHE':'F', 'TRP':'W', 'TYR':'Y',
			   'GLY':'G', 'ALA':'A', 'ILE':'I', 'LEU':'L', 'CYS':'C',
			   'MET':'M', 'THR':'T', 'VAL':'V', 'PRO':'P',
			   'HID':'H', 'HIE':'H', 'HIP':'H',
			   'ASH':'D', 'GLH':'E',
			   'LYN':'K', 'ARN':'R',
			   'HOH':'U', 'CL': 'J' }

ADbatchRE = r'FAHV_(x?)(.+)_([0-9]+)_processed.tgz'
ADbatchREPat = re.compile(ADbatchRE)

InterTypes = ('hba', 'hbd', 'mtl','ppi','tpi','vdw')
InterRE = r'(.*):.+():([A-Za-z]+[0-9]*)~~(.*):(.+):(.+)'
InterREPat = re.compile(InterRE)

# vdw are different
# :CYS65:SG
# InterVDWRE = r'(.*):([A-Z]+[0-9]+):([A-Z]+)'
InterVDWRE = r'(.*):(.+):(.+)'
InterVDWREPat = re.compile(InterVDWRE)

# ppi/tpi are different
# B:HIS114~~(-4.295,-13.390,-20.427:-3.408,-10.290,-18.368)
# cf. piStackingAndRingDetection.findLigRecPiStack()
# pstack.append([res, rec_centroid, lig_centroid])
InterPiRE = r'(.*):(.+)~~\(([-.,0-9]+):([-.,0-9]+)\)'
InterPiREPat = re.compile(InterPiRE)

def reducePlusInterDict(ligData):
	'''create sparse version of interaction dict, with only ligands atom type (not its index)
	   the syntax for the InterPattern seems a PDB convention?
	'''
	
	rinterDict = {}
	for itype in InterTypes:
		if len(ligData[itype]) > 0:
			rinterDict[itype] = []
			if itype=='ppi' or itype=='tpi':
				for inter in ligData[itype]:
					# d:<0>:O3~~B:ARG57:N
					m = InterPiREPat.match(inter)
					try:
						(rchain,raa,rcenter,ligcenter) = m.groups()
						rinterDict[itype].append( (rchain,raa,rcenter,ligcenter) )
					except:
						print 'reducePlusInterDict: bad ppi/tpi string?!',inter
			else:
				# 151110: new processing now includes vdw ligand atoms
				for inter in ligData[itype]:
					# d:<0>:O3~~B:ARG57:N
					m = InterREPat.match(inter)
					try:
						# (ligIdx,liname,rchain,raa,ratom) = m.groups()
						(lchain,lres,latom,rchain,raa,ratom) = m.groups()
						rinterDict[itype].append( (rchain,raa,ratom,latom) )
					except:
						print 'reducePlusInterDict: bad %s string?! %s' % (itype,inter)
			
	return rinterDict

#############################################################################
# from rabbit
# Author: Stefano FORLI
#
# Copyright: Stefano Forli, TSRI 2011
#
# v.0.4
#############################################################################

def checkVSresult_ADV(lines):
	"parses and validates a PDBQT+ file (lines)"
	if (lines[0].startswith("USER    ADVS_Vina_result>")):
		return True
	else:
		return False

def setKwMode_ADV(mode = "1"):
	# USER ADVina_pose1> -10.200, -0.309
	# mode = "1", "2", "..."
	if mode == "any":
		kw = "ADVina_pose."
	else:
		kw = "ADVina_pose"+mode
	return kw

def getResultCount_ADV(lines):
	return int(lines[4].split("ADVina_results>")[1].strip())
	
def getRawEnergyData_ADV(lines, mode = "1"):
	# mode = "1", "2", "...", "any" (POSES)
	# format: { "e" : [ float(e)] , "leff" : [float(leff)] }
	kw = setKwMode_ADV(mode)+">"
	result = { "e" : [],
				"leff" : [] }
	for l in lines:
		if re.search(kw, l):
			l = l.split(">", 1)[1]
			e, leff = l.split(",")
			result['e'].append(float(e))
			result["leff"].append(float(leff))
			if not mode == "any":
				break
	return result

def getLigInteractions_ADV(lines, mode = "1"):
	kw = setKwMode_ADV(mode) #
	kw += "_" # TODO test this!
	interactions = {"vdw" : [], # these keys must match the
					"hba" : [], # tags used in writing the PDBQT+
					"hbd" : [],
					"ppi" : [],
					"tpi" : [],
					"mtl" : [],
					}
	for l in lines:
		if re.search(kw, l):
			for itype in interactions.keys():
				if (itype == "ppi") or (itype == "tpi"):
					sep = ";"
				else:
					sep = ","
				if itype in l:
					l = l.split(itype+">", 1)[1]
					l = l.split(sep)
					for i in l:
						interactions[itype].append(i.strip())
	return interactions

def getLigSource_ADV(lines):
	srcPrefix = 'USER    ADVina_pose'
	src = ''
	## 160531: grabs FIRST pose number in PDBQT?!
	for l in lines:
		if l.startswith(srcPrefix):
			gtpos = l.find('>')
			src = l[gtpos+1:].strip()
			return src
	return src

def getGenericData_ADV(lines):
	receptNamePrefix = 'USER    ADVina_rec> '
	nresultPrefix = 'USER    ADVina_results> '
	rname = ''
	nresult = 0
	for l in lines:
		if l.startswith(receptNamePrefix):
			rname = l[len(receptNamePrefix):].strip()
		if l.find(nresultPrefix) != -1:
			nresult = int(l[len(nresultPrefix):].strip())

	return (rname,nresult)

#############################################################################
# eo rabbit
	
def parseADPDBQT_ADV(f):
	'''uses raccoon's routines to extract 
		{'recept': '.../structure/x3NF6_B_IN_LEDGF.pdbqt',
		 'e': -6.9,
		 'leff': -0.216,
		 'nresult': 1,
		 'src': 's> 1',
		 'hba': [('A', 'GLU170', 'N', 'O33'), ...],
		 'vdw': [('A', 'GLU170', 'HN', 'C25'), ... ]
		 }
	'''
	ligand = getLines(f)
	if not checkVSresult_ADV(ligand):
		return None

	ligData = {}
	rname,nresult = getGenericData_ADV(ligand)
	
	if rname=='':
		print 'parseADPDBQT_ADV: missing receptor name?!',f
		ligData['recept'] = ''
	else:
		ligData['recept'] = rname
	
	if nresult==0:
		print 'parseADPDBQT_ADV: missing results?!',f
		ligData['nresult'] = 0
		return None
	else:
		ligData['nresult'] = nresult
		
	## ADVina has one
	src = getLigSource_ADV(ligand)
	if src=='':
		print 'parseADPDBQT_ADV: missing src?!',f
		ligData['nresult'] = 0
	else:
		ligData['src'] = src
	
	ligdataRaw = getRawEnergyData_ADV(ligand)
	# dict: {'c_pc': [53.45], 'e': [-6.19], 'leff': [-0.413], 'c_size': [93]}

	for k,v in ligdataRaw.items():
		ligData[k] = v[0]
		
	liginteract = getLigInteractions_ADV(ligand)
	reducedInterDict = reducePlusInterDict(liginteract)
	ligData.update(reducedInterDict)
	
	# do additional pdbqt access via OB
	# to convert tpi, ppi ligandCenter to closest atom
	if 'ppi' in ligData or 'tpi' in ligData:
	
		allMol = pybel.readfile('pdbqt', f)
		# print 'len(allMol)',len(allMol)
		ligMol = allMol.next() # ASSUME only one mol PDBQT
		obmol = ligMol.OBMol
		atomCoord = {}
		for res in ob.OBResidueIter(obmol):
			for obatom in ob.OBResidueAtomIter(res):
				pbatom = pybel.Atom(obatom)
				idx = pbatom.idx
				laName = res.GetAtomID(obatom).strip()
				laIdx = idx
				laCoord = pbatom.coords
				laType = pbatom.type
				atomCoord[idx] = [laName,laType,laCoord]
		
		for itype in ['tpi' ,'ppi']:
			if itype not in ligData:
				continue
			
			newInterList = []
			for inter in ligData[itype]:
				#NB: ratom is really receptorCenter for tpi/ppi
				(rchain,raa,rcenter,lcenter) = inter
				minLA = None
				minDist = 1.0e6
				laTuple = tuple(eval(lcenter)) # NB: latomFull is a string
				obLAnp = np.array(laTuple)
				for laIdx in atomCoord:
					laName,laType,laCoord = atomCoord[laIdx]
					latype = FAAHA.getAtomType(laName)
					# 151109: only consider ligatoms *OTHER THAN* H, C
					if latype in config.VDWExcludedLigAtoms:
						continue
					  
					rlifLAnp = np.array(laCoord)
					l2norm = np.linalg.norm(obLAnp - rlifLAnp)
					# print laIdx,atName,l2norm
					if l2norm < minDist:
						minLA = laIdx
						minDist = l2norm
						
				ratom = ''
				latomFull = atomCoord[minLA][0]
				newInterList.append( (rchain,raa,ratom,latomFull) )		  
			ligData[itype] = newInterList
				   
	return ligData
		

# simplified, here, means that this should work for all batches after ??? Exp.# ???
# because naming was simplified (by dns)
# ADVsimple = r'fahv.x([a-zA-Z0-9]+)(_[A-Z0-9]*)*_(ZINC[0-9]+)(_[0-9]*)*_([0-9]+)_out_Vina_VS.pdbqt'
# fahv.x3KF0_prASw0c0_ZINC00147966_544588015_out_Vina_VS.pdbqt

# 160110: ligand capturing part of qualified receptor?
# fahv.x3KF0_prASw0c0_ZINC00147966_544588015_out_Vina_VS.pdbqt
# ADVsimple = r'fahv.x([a-zA-Z0-9]+)_(.+)_([0-9]+)_out_Vina_VS.pdbqt'
ADVsimple1 = r'fahv.x([a-zA-Z0-9_]+)_(.+)_([0-9]+)_out_Vina_VS.pdbqt'
ADVsimplePat1 = re.compile(ADVsimple1, re.IGNORECASE)

# 160216: consuming SForli's rerun of DUDE
# ZINC63694490_out_Vina_VS.pdbqt
# activated by config.procVersion = 0.2
ADVsimple2 = r'([A-Z0-9_]+)_out_Vina_VS.pdbqt'
ADVsimplePat2 = re.compile(ADVsimple2, re.IGNORECASE)

def visit_ADV_tgz(tgzPath,exptname,recon,verbose):

	dataTbl = {}
	# 2do: Py2.7 allows WITH context management! TODO
# with tarfile.open(tgzPath) as subTar:
# with tarfile.open(subTar) as dataDir:

	tmpDir = tempfile.mkdtemp()
	allTar = tarfile.open(tgzPath)
	allTar.extractall(tmpDir)

	# ASSUME: _VS "bar" style processed file names for ADV
	# fahv.x4I7G_RT_NNRTIadj_wNNRTI_ZINC58421065_1_649284996_out_Vina_VS.pdbqt
	# Exp79/Results_x3kf0A/FAHV_x3kf0A_0124403_processed.tgz example:
	#  FAHV_x3kf0A_0124403_processed/ # (untarred dir)
	#	fahv.x3kf0A_ZINC01569654_1113915765_out_Vina_VS.pdbqt

	procList = glob.glob(tmpDir+'/FAHV*/fahv.*_out_Vina_VS.pdbqt')
	if verbose:
		print 'visit_ADV_tgz: NTGZ=',len(procList)
		
	ndup =0
	ninvFile=0
	for isd,procPath in enumerate(procList):
		
		procBits = os.path.split(procPath)

		# 
		dirBits = procBits[0].split('/')
		currReceptor = dirBits[-2]
		
		# batchNo = int(dirBits[-1].split('_')[1])

		procf = procBits[1]
		
		# lbpos = procf.find('_')
		# rbpos = procf.rfind('_')
		# assert (lbpos != -1 and rbpos != -1 and rbpos > lbpos), 'visit_ADV_tgz: bad procf?! %s' % (procf)
		# ligand = procf[lbpos+1:rbpos]
		
		if config.procVersion == '0.2':
			mpath = ADVsimplePat2.match(procf)
			ligand = mpath.groups()[0]
			receptor = currReceptor
		else:
			mpath = ADVsimplePat1.match(procf)
			(receptor,ligand,workNo) = mpath.groups()
				
		dk = (exptname,receptor,ligand)
		
		if dk in dataTbl:
			print 'visit_ADV_tgz: dup dataKey?!',dk
			ndup += 1
			continue

		###-------
		ligData = parseADPDBQT_ADV(procPath)
		###-------
		
		if not(ligData):
			print 'visit_ADV_tgz: invalid ADV file?!',procf, tgzPath
			ninvFile += 1
			continue

		dataTbl[dk] = ligData
		
		if verbose and isd % 100 == 0:
			print 'visit_ADV_tgz: nproc=', isd
						
	shutil.rmtree(tmpDir)
	print 'visit_ADV_tgz: done. NLig=%d NDup=%d NInvalidFile=%d' % \
		(len(dataTbl),ndup,ninvFile)
	
	return dataTbl

def visit_ADV(procPath,exptname,recon,verbose):
	'''assume path already has processed files (aot/ tgzPath), as produced by process
	'''

	dataTbl = {}
	# 2do: Py2.7 allows WITH context management! TODO
# with tarfile.open(tgzPath) as subTar:
# with tarfile.open(subTar) as dataDir:

	# ASSUME: _VS "bar" style processed file names for ADV
	# fahv.x4I7G_RT_NNRTIadj_wNNRTI_ZINC58421065_1_649284996_out_Vina_VS.pdbqt
	# Exp79/Results_x3kf0A/FAHV_x3kf0A_0124403_processed.tgz example:
	#  FAHV_x3kf0A_0124403_processed/ # (untarred dir)
	#	fahv.x3kf0A_ZINC01569654_1113915765_out_Vina_VS.pdbqt

	procList = glob.glob(procPath+'/*_out_Vina_VS.pdbqt')
	if verbose:
		print 'visit_ADV: NProcFiles=',len(procList)
		
	ndup =0
	ninvFile=0
	for isd,procPath in enumerate(procList):
		# fahv.x3kf0A_ZINC00145439_2057149382_out_Vina_VS.pdbqt
		procBits = os.path.split(procPath)

		# 
		dirBits = procBits[0].split('/')
		currReceptor = dirBits[-2]
		
		# batchNo = int(dirBits[-1].split('_')[1])

		procf = procBits[1]

		# 160110: ligand capturing part of qualified receptor?
		# fahv.x3KF0_prASw0c0_ZINC00147966_544588015_out_Vina_VS.pdbqt

		# fahv.x3kf0A_ZINC00145439_2057149382_out_Vina_VS.pdbqt
		
		# ADVsimple = r'fahv.x([a-zA-Z0-9]+)_(.+)_([0-9]+)_out_Vina_VS.pdbqt'
		
		# lbpos = procf.find('_')
		# rbpos = procf.rfind('_')
		# assert (lbpos != -1 and rbpos != -1 and rbpos > lbpos), 'visit_ADV_tgz: bad procf?! %s' % (procf)
		# ligand = procf[lbpos+1:rbpos]
		
		if config.procVersion == '0.2':
			mpath = ADVsimplePat2.match(procf)
			ligand = mpath.groups()[0]
			receptor = currReceptor
		else:
			mpath = ADVsimplePat1.match(procf)
			(receptor,ligand,workNo) = mpath.groups()
		
		dk = (exptname,receptor,ligand)
		
		if dk in dataTbl:
			print 'visit_ADV: dup dataKey isd=%d %s?!' % (isd,dk)
			ndup += 1
			continue

		###-------
		ligData = parseADPDBQT_ADV(procPath)
		###-------
		
		if not(ligData):
			print 'visit_ADV: invalid ADV file?!',procf, procPath
			ninvFile += 1
			continue

		dataTbl[dk] = ligData

# 		if verbose and isd % 1000 == 0:
# 			print 'visit_ADV: nproc=', isd

	print 'visit_ADV: procPath=%s done. NLig=%d NDup=%d NInvalidFile=%d' % \
		(procPath,len(dataTbl),ndup,ninvFile)
	
	return dataTbl

def rptData_ADV(dataTbl,summf,interf,exptname,batchNo):
	'''V1: produce condensed JSON inter file
			[ [Expt,BatchNo,Recept,Lig, [IType,[InterEnum] ] ] ]
			ala ["Exp96", 197339, "x3ZSW_B_IN_Y3", "Y3_ZINC00626007", [[0, [["B", "R199", "NH2", "O1"], ["B", "K188", "NZ", "O2"]]], [5, [["B", "G82", "CA"], ..., ["B", "I141", "O"]]]]]
	'''
	
	summs = open(summf,'w')
	summs.write('Expt,Batch,Recept,Ligand,E,Eff,Nvdw,Ninter\n')
	
	allInter = []
	for dk in dataTbl:
		(exptname,receptor,lig) = dk
		ligData = dataTbl[dk]
		ninter = 0
		nvdw = 0
		for itype in InterTypes:
			if itype in ligData:
				if itype=='vdw':
					nvdw = len(ligData['vdw'])
				else:
					ninter += len(ligData[itype])

		summs.write('%s,%d,%s,%s,%s,%s,%d,%d\n' % \
				   (exptname,batchNo,ligData['recept'],lig,\
					ligData['e'],ligData['leff'],nvdw,ninter))
				

		interInfo = [exptname, batchNo, ligData['recept'], lig]
		interList = []
		
		for itype in InterTypes:
			if itype in ligData:
				# additional compression by combining all inter of same type
				itlist = []
				# convert itype string to its index in InterTypes
				itypeIdx = InterTypes.index(itype)

				for interTuple in ligData[itype]:
					inter = list(interTuple)
					
					# convert 3-letter receptor AA to single char
					raa = inter[1]
					raalet3 = raa[:3]
					pos = raa[3:]
					if raalet3 in AADict:
						raaLet = AADict[raalet3]
					else:
						raaLet = 'X'
					inter[1] = raaLet+pos
					
					itlist.append(inter)
					
				interList.append( [itypeIdx,itlist] )
					
		interInfo.append(interList)
		
		allInter.append(interInfo)
		
	summs.close()
				  
	inters = open(interf,'w')
	json.dump(allInter,inters)
	inters.close()


def visitRpt_ADV_tgz(tgzPath,outPath,recon,batchTbl,tocs,exptname,batchNo,verbose):
	''' get info, place in this table
		info is extracted from enhanced pdbqt files inside tarball
	'''	

	dataTbl = visit_ADV_tgz(tgzPath,exptname,recon,verbose)
	if recon:
		print 'visitRpt_ADV_tgz: Recon-only; no reporting'
		return len(dataTbl)
		
	summf  = outPath +'/summ/ADV_summ_%07d.csv' % (batchNo)
	interf = outPath +'/inter/ADV_inter_%07d.json' % (batchNo)
	rptData_ADV(dataTbl,summf,interf,exptname,batchNo)
	
	#tocStr = '%s,%05d,%d,%s' % (exptname,ntgz,len(dataTbl),tgzPath)
	tocStr = '%s,%d,%d,%s' % (exptname,batchNo,len(dataTbl),tgzPath)
	if verbose:
		print 'visitRpt_ADV_tgz toc:',tocStr
			
	tocs.write(tocStr+'\n')
	# fun2watch! toc
	tocs.flush(); os.fsync(tocs.fileno())
	
	return len(dataTbl)

def visitRpt_ADV(tgzPath,outPath,recon,batchTbl,tocs,exptname,batchNo,verbose):
	''' get info, place in this table
		info is extracted from enhanced pdbqt files inside tarball
		assume path already has processed files (aot/ tgzPath), as produced by process
	'''	

	dataTbl = visit_ADV(tgzPath,exptname,recon,verbose)
	if recon:
		print 'visitRpt_ADV: Recon-only; no reporting'
		return len(dataTbl)
		
	summf  = outPath +'summ/ADV_summ_%07d.csv' % (batchNo)
	interf = outPath +'inter/ADV_inter_%07d.json' % (batchNo)
	rptData_ADV(dataTbl,summf,interf,exptname,batchNo)
	
	#tocStr = '%s,%05d,%d,%s' % (exptname,ntgz,len(dataTbl),tgzPath)
	tocStr = '%s,%d,%d,%s' % (exptname,batchNo,len(dataTbl),tgzPath)
	if verbose:
		print 'visitRpt_ADV toc:',tocStr
			
	tocs.write(tocStr+'\n')
	# fun2watch! toc
	tocs.flush(); os.fsync(tocs.fileno())
	
	return len(dataTbl)

def mglTop_visit_ADV(tgzProc=True,exptList=None,recon=False,verbose=False):
	'recon stops after opening, parsing one file in first tgz'
	
	if exptList:
		print 'mglTop_visit_ADV: Explicit experiment list %s' % (str(exptList))
	else:
		print 'mglTop_visit_ADV: Full crawl of %s' % (ProcDir)
		dirbits = ProcDir.split('/')
		# ASSUME ADV_topDir = .../ exptName / receptor /
		#			ala       .../DUDE-HIVPR/x1XL2/
		pexpt = dirbits[-3]
		
		# exptList = [os.path.split(exptPath)[1] for exptPath in glob.glob(crawlPat) ]
		# 160218: Empty string causes exptPath = ADV_topDIr, outPath = outDir
		exptList = [pexpt]
		
	print 'mglTop_visit_ADV: NExperiments=',len(exptList)
	
	if verbose:
		print 'mglTop_visit_ADV: **Verbose output'

	if recon:
		print 'mglTop_visit_ADV: **Reconnaissance sniffing only!'

	totParse = 0
	for ie, exptname in enumerate(exptList):

		startTime = datetime.datetime.now()
		print 'mglTop_visit_ADV: %s starting %s' % (exptname,startTime.strftime("%y%m%d %H:%M:%S"))

		# only create exptname subdir if exptname does NOT come from project experiment, above
		if ProcDir.find(exptname) == -1:
			exptPath = ProcDir+exptname+'/'
			outPath = CrawlDir+exptname+'/'
		else:
			exptPath = ProcDir
			outPath = CrawlDir
		
		if verbose:
			print ' *',ie,exptPath
			
		exptReceptList = glob.glob(exptPath+'*')
		exptReceptList.sort()
		for receptPath in exptReceptList:
			recept = receptPath.split('/')[-1]
				
			exptReceptPath = exptPath + recept + '/'
			outReceptPath = outPath + recept + '/'

			if verbose:
				print ' **',ie,recept,exptReceptPath

			if not os.path.isdir(outReceptPath):
				print 'mglTop_visit_ADV: creating ExptOutput directory',outReceptPath
				os.makedirs(outReceptPath)
					
			tocf = outReceptPath+'ADV_toc.csv'
			tocs = open(tocf,'w')
			#tocs.write('NTGZ,Data,Path\n')
			tocs.write('Experiment, Batch, Data, Path\n')
	
			summDir = outReceptPath+'summ/'
			if not os.path.isdir(summDir):
				print 'mglTop_visit_ADV: creating summDir ',summDir
				os.makedirs(summDir)
	
			interDir = outReceptPath+'inter/'
			if not os.path.isdir(interDir):
				print 'mglTop_visit_ADV: creating interDir ',interDir
				os.makedirs(interDir)
	
			batchTbl = {}
	
			if tgzProc:
				exptSubList = glob.glob(exptReceptPath+'batch_*.tgz')
			else:
				exptSubList = glob.glob(exptReceptPath+'batch_*')   
	
			print 'mglTop_visit_ADV: Recept=%s NBatches=%d' % (recept,len(exptSubList))
			#NB: glob() returns matches in arbitrary order
			exptSubList.sort()
			for ise, exptSubPath in enumerate(exptSubList):
				# exptSubPath =  /Data/sharedData/coevol-HIV/WCG/process2/Exp120/Results_x3KF0_prASw0c0/batch_0550340
				if verbose:
					print ' ***',ie,ise,recept,exptSubPath
	
				if tgzProc:
					tgzList = glob.glob(exptSubPath+'/*.tgz') 
					print 'mglTop_visit_ADV: NTGZ=',len(tgzList)
					tgzList.sort()
					for jt,tgzPath in enumerate(tgzList):
						if recon and jt > 0:
							print 'mglTop_visit_ADV: Recon-only; break'
							break
											
							tgznow = os.path.split(tgzPath)[1]
							mpath = ADbatchREPat.match(tgznow)
							if mpath:
								(x, vinarec, vinabatch) = mpath.groups()
							else:
								print 'mglTop_visit_ADV: bad match?!',tgznow
								continue
							batchNo = int(vinabatch)
							if verbose:
								print 'mglTop_visit_ADV: Attempting to analyze',tgznow 
		
							nparse = visitRpt_ADV_tgz(tgzPath,outReceptPath,recon,batchTbl,tocs,exptname,batchNo,verbose)
				else:
					# exptSubPath = /Data/sharedData/coevol-HIV/WCG/process2/Exp120/batch_0550340
					
					bpos = exptSubPath.rfind('_')
					bnos = exptSubPath[bpos+1:]
					batchNo = int(bnos)
					
					nparse = visitRpt_ADV(exptSubPath,outReceptPath,recon,batchTbl,tocs,exptname,batchNo,verbose)
						#
				totParse += nparse
					
			endTime = datetime.datetime.now()
			elapTime = endTime-startTime
			print 'mglTop_visit_ADV: Expt=%s Recept=%s done. TotParse=%d NSec=%s' % (exptname,recept,totParse,elapTime.total_seconds())
				
			tocs.close() # for each experiment directory
		
	print 'mglTop_visit_ADV: TotParse=',totParse

def Local_Top_visit_ADV2(ADV_topDir,outdir,exptName,receptor,procPattern = '*_out_Vina_VS.pdbqt',tgzProc=False,verbose=False):
	'crawl local (non-FAAH) processed files'

	tocf = outdir+'ADV_toc.csv'
	tocs = open(tocf,'w')
	#tocs.write('NTGZ,Data,Path\n')
	tocs.write('Experiment, Batch, Data, Path\n')

	summDir = outdir+'summ/'
	if not os.path.isdir(summDir):
		print 'mglTop_visit_ADV: creating summDir ',summDir
		os.makedirs(summDir)

	interDir = outdir+'inter/'
	if not os.path.isdir(interDir):
		print 'mglTop_visit_ADV: creating interDir ',interDir
		os.makedirs(interDir)

	batchTbl = {}

	if tgzProc:
		batchList = glob.glob(ADV_topDir+'batch_*.tgz')
	else:
		batchList = glob.glob(ADV_topDir+'batch_*')   
		
	print 'Local_Top_visit_ADV2: NBatches=',len(batchList)
	#NB: glob() returns matches in arbitrary order
	batchList.sort()
	
	totParse = 0
	recon = False
	
	for ib, batchPath in enumerate(batchList):
		
		bpos = batchPath.rfind('_')
		bnos = batchPath[bpos+1:]
		batchNo = int(bnos)
		
		nparse = visitRpt_ADV(batchPath,recon,batchTbl,outdir,tocs,exptName,batchNo,verbose)
				#
		totParse += nparse

# 		procList = pathToList(ADV_topDir,recursive=True, pattern=procPattern)
# 		
# 		totParse = 0
# 		dataTbl = {}
# 		print 'Local_Top_visit_ADV2: NProc=',len(procList)
# 		for isd,procPath in enumerate(procList):
# 	
# 			procf = os.path.split(procPath)[1]
# 			fbits = procf.split('_')
# 			ligand = fbits[0]
# 					
# 			###-------
# 			ligData = parseADPDBQT_ADV(procPath)
# 			###-------
# 	
# 			if not(ligData):
# 				print 'Local_Top_visit_ADV2: invalid ADV file?!',procPath
# 				continue
# 	
# 			dk = (exptName,receptor,ligand)
# 			
# 			if dk in dataTbl:
# 				print 'Local_Top_visit_ADV2: dup dataKey?!',dk
# 				continue
# 	
# 			dataTbl[dk] = ligData
# 								   
# 			totParse += 1
#
# 	summf  = outdir+'/ADV_summ.csv'
# 	interf = outdir+'/ADV_inter.json' 
# 	rptLocalData_ADV(dataTbl,summf,interf)

	print 'Local_Top_visit_ADV2: TotParse=',totParse
			
# 141028 rikPR253
# .../processed/rikPR253/x1A8K_PRAS/PDB0Q4_out_Vina_VS.pdbqt
Local_PR253_RE = r'(.+)_out_Vina_VS.pdbqt'
Local_PR253_Pat = re.compile(Local_PR253_RE)

def Local_PR253_Top_visit_ADV(ADV_topDir,outdir,verbose=False):
	'crawl local (non-FAAH) processed files'
	
	# .../processed/rikPR253/x1A8K_PRAS/PDB0Q4_out_Vina_VS.pdbqt

	totParse = 0
	dataTbl = {}
	pdbDirList = glob.glob(ADV_topDir+'/x*_PRAS') 
	print 'Local_PR253_Top_visit_ADV: NSubExperiments=',len(pdbDirList)
	for ise, exptSubPath in enumerate(pdbDirList):
		
		pathBits = exptSubPath.split('/')
		receptor = exptSubPath.split('/')[-1]
		assert receptor.endswith('_PRAS')
		assert receptor[0]=='x'
		receptor = receptor[1:]
		receptor = receptor[:-5]
		
		procList = glob.glob(exptSubPath+'/*_out_Vina_VS.pdbqt')
		for isd,procPath in enumerate(procList):

			procf = os.path.split(procPath)[1]
			mpath = Local_PR253_Pat.match(procf)
			if mpath:
				ligand = mpath.groups()[0]
			else:
				print 'Local_PR253_Top_visit_ADV: bad match?!',procPath
				continue
 
			###-------
			ligData = parseADPDBQT_ADV(procPath)
			###-------

			if not(ligData):
				print 'Local_PR253_Top_visit_ADV: invalid ADV file?!',procPath
				continue
	
			dk = (receptor,ligand)
			
			if dk in dataTbl:
				print 'Local_PR253_Top_visit_ADV: dup dataKey?!',dk
				continue
	
			dataTbl[dk] = ligData
								   
			totParse += 1
			
	print 'Local_PR253_Top_visit_ADV: TotParse=',totParse
			
	summf  = outdir+'/ADV_summ.csv'
	interf = outdir+'/ADV_inter.json' 
	rptLocalData_ADV(dataTbl,summf,interf)

# 141016 rikPR19 dockings: 19 inhib X 
# .../wcg/processed/rikPR19/processed/Results_x1A30-AS/PDB017_x1A30-AS_vina-out_Vina_VS.pdbqt
#								  ...Results_x1FB7-1w0n0-AAsh25/PDB017_x1FB7-1w0n0-AAsh25_vina-out_Vina_VS.pdbqt
#
LocalRE_PR19 = r'(.+)_(.+)_vina-out_Vina_VS.pdbqt'
LocalRE_PR19_Pat = re.compile(LocalRE_PR19)

def PR19_Top_visit_ADV(ADV_topDir,outdir,verbose=False):
	'crawl local (non-FAAH) processed files'
	
	# 'Results_%s/%s_%s_vina-out_VS.pdbqt' % (recept,lig,recept)

	totParse = 0
	dataTbl = {}
	exptSubList = glob.glob(ADV_topDir+'/Results_*') 
	print 'PR19_Top_visit_ADV: NSubExperiments=',len(exptSubList)
	for ise, exptSubPath in enumerate(exptSubList):
		if verbose:
			print ' **',ise,exptSubPath
			
		procList = glob.glob(exptSubPath+'/*_vina-out_Vina_VS.pdbqt')
		for isd,procPath in enumerate(procList):
			# fahv.x3kf0A_ZINC00145439_2057149382_out_Vina_VS.pdbqt

			procf = os.path.split(procPath)[1]
			mpath = LocalRE_PR19_Pat.match(procf)
			if mpath:
				(ligand,receptor) = mpath.groups()
			else:
				print 'PR19_Top_visit_ADV: bad match?!',procPath
				continue
 
			###-------
			ligData = parseADPDBQT_ADV(procPath)
			###-------

			if not(ligData):
				print 'PR19_Top_visit_ADV: invalid ADV file?!',procPath
				continue
	
			dk = (receptor,ligand)
			
			if dk in dataTbl:
				print 'PR19_Top_visit_ADV: dup dataKey?!',dk
				continue
	
			dataTbl[dk] = ligData
								   
			totParse += 1
			
	print 'PR19_Top_visit_ADV: TotParse=',totParse
			
	summf  = outdir+'/ADV_summ.csv'
	interf = outdir+'/ADV_inter.json' 
	rptLocalData_ADV(dataTbl,summf,interf)
			  
def rptLocalData_ADV(dataTbl,summf,interf):
	'''V1: produce condensed JSON inter file
			[ [Expt,BatchNo,Recept,Lig, [IType,[InterEnum] ] ] ]
			ala ["Exp96", 197339, "x3ZSW_B_IN_Y3", "Y3_ZINC00626007", [[0, [["B", "R199", "NH2", "O1"], ["B", "K188", "NZ", "O2"]]], [5, [["B", "G82", "CA"], ..., ["B", "I141", "O"]]]]]
	'''
	
	summs = open(summf,'w')
	summs.write('Recept,Ligand,E,Eff,Nvdw,Ninter\n')
	
	allInter = []
	for dk in dataTbl:
		(exptname,receptor,lig) = dk
		ligData = dataTbl[dk]
		ninter = 0
		nvdw = 0
		for itype in InterTypes:
			if itype in ligData:
				if itype=='vdw':
					nvdw = len(ligData['vdw'])
				else:
					ninter += len(ligData[itype])

		summs.write('%s,%s,%s,%s,%d,%d\n' % \
				   (receptor,lig,\
					ligData['e'],ligData['leff'],nvdw,ninter))
				

		interInfo = [ receptor, lig]
		interList = []
		
		for itype in InterTypes:
			if itype in ligData:
				# additional compression by combining all inter of same type
				itlist = []
				# convert itype string to its index in InterTypes
				itypeIdx = InterTypes.index(itype)

				for interTuple in ligData[itype]:
					inter = list(interTuple)
					
					# convert 3-letter receptor AA to single char
					raa = inter[1]
					raalet3 = raa[:3]
					pos = raa[3:]
					if raalet3 in AADict:
						raaLet = AADict[raalet3]
					else:
						raaLet = 'X'
					inter[1] = raaLet+pos
					
					itlist.append(inter)
					
				interList.append( [itypeIdx,itlist] )
					
		interInfo.append(interList)
		
		allInter.append(interInfo)
		
	summs.close()
				  
	inters = open(interf,'w')
	json.dump(allInter,inters)
	inters.close()

def Local_focusedLib_Top_visit_ADV(fldir,outdir,verbose=False):
	'crawl Stefanos "focused library" of ILINI compounds'
	
	# .../processed/focusedLib/lib_v2.1__R1/chlorine/0/libR1_clorinated_5/libR1_clorinated_5_out_Vina_VS.pdbqt
	
	print 'Local_focusedLib_Top_visit_ADV: NFiles=%d (_VS=%d .VS=%d)' % (len(allFiles),len(fileList1),len(fileList2))
	
	totParse = 0
	dataTbl = {}
	exptname = 'focusLib'
	receptor = 'CCDKF115'
	fpathFile = outdir+'filePaths.csv'
	outs = open(fpathFile,'w')
	outs.write('Lib,Ligand,Path\n')
	allPaths = {}
	for fi, path in enumerate(allFiles):
		
		pathBits = path.split('/')
		assert pathBits[6] == 'focusedLib', 'Odd path?!'
		lib = pathBits[7]
 
		###-------
		ligData = parseADPDBQT_ADV(path)
		###-------

		if not(ligData):
			print 'Local_focusedLib_Top_visit_ADV: invalid ADV file?!',path
			continue
	
		fname = pathBits[-1]
		if fname.find('_VS') != -1:
			fbits = fname.split('_')
			try:
				outPos = fbits.index('out')
			except:
				print 'huh'
			if fbits[0].startswith('lib'):
				# libR1_clorinated_5_out_Vina_VS.pdbqt
				ligand = '_'.join(fbits[1:outPos])
			else:
				# ligand435_mini_ph_out_Vina_VS.pdbqt
				ligand = '_'.join(fbits[:outPos])			
			
		elif fname.find('.VS') != -1:
			# lib01_27_minimized_ph.VS.pdbqt
			ppos = fname.find('.')
			fbits = fname[:ppos].split('_')
			# NB: all of these lib01 files have lib prefix
			ligand = '_'.join(fbits[1:])
		
		assert ligand not in allPaths, 'Local_focusedLib_Top_visit_ADV: dup ligand?! %s %s' % (lib,ligand)
		allPaths[ligand] = path
		relPath = path.replace(ProcDir,'')
		outs.write('%s,%s,%s\n' % (lib,ligand,relPath))
		
		dk = (exptname,receptor,ligand)
		
		if dk in dataTbl:
			print 'Local_focusedLib_Top_visit_ADV: dup dataKey?!',dk
			continue

		dataTbl[dk] = ligData
							   
		totParse += 1
	
	outs.close()		
	print 'Local_focusedLib_Top_visit_ADV: TotParse=',totParse
	
	# NB: mimic experiment, batch numbering scheme
	summf  = outdir+ ('Exp1/summ/ADV_summ_%07d.csv' % (1))
	interf = outdir+ ('Exp1/inter/ADV_inter_%07d.json' % (1))
	rptLocalData_ADV(dataTbl,summf,interf)
	

def Local_Top_visit_ADV(fldir,outdir,exptname,receptName,fileList,verbose=False):
	'crawl local SAMPL4 libraries of processed results'
	
	# .../processed/focusedLib/lib_v2.1__R1/chlorine/0/libR1_clorinated_5/libR1_clorinated_5_out_Vina_VS.pdbqt
		
	totParse = 0
	dataTbl = {}
	ndup = 0
	fpathFile = outdir+'filePaths.csv'
	outs = open(fpathFile,'w')
	outs.write('Expt,Lib,Recept,Ligand,Path\n')
	allPaths = {}
	for fi, path in enumerate(fileList):
		
		pathBits = path.split('/')
 
		###-------
		ligData = parseADPDBQT_ADV(path)
		###-------

		if not(ligData):
			print 'Local_Top_visit_ADV: invalid ADV file?!',path
			continue
	
		if RunName == 'focusedLib':
			fname = pathBits[-1]
			lib = pathBits[7]
			if fname.find('_VS') != -1:
				fbits = fname.split('_')
				try:
					outPos = fbits.index('out')
				except:
					print 'huh'
				if fbits[0].startswith('lib'):
					# libR1_clorinated_5_out_Vina_VS.pdbqt
					ligand = '_'.join(fbits[1:outPos])
				else:
					# ligand435_mini_ph_out_Vina_VS.pdbqt
					ligand = '_'.join(fbits[:outPos])			
				
			elif fname.find('.VS') != -1:
				# lib01_27_minimized_ph.VS.pdbqt
				ppos = fname.find('.')
				fbits = fname[:ppos].split('_')
				# NB: all of these lib01 files have lib prefix
				ligand = '_'.join(fbits[1:])
		elif RunName == 'SAMPL4':
			# pathBits = ['', 'Data', 'sharedData', 'coevol-HIV', 'WCG', 'processed', 'SAMPL4', 'LEDGF', '', 's3NF8_A', 'AVX101118_0', 'AVX101118_0_VINoutput_Vina_VS.pdbqt']
			ligand = pathBits[-2]
			lib = RunName
			
		else:
			# DUDE-151110
			fname = pathBits[-1]
			# HACK: ASSUME ligand ends with digits
			for i in range(len(fname)-1,0,-1):
				if fname[i].isdigit():
					eoLig = i 
					break
			ligand = fname[:eoLig+1]
			lib = RunName 
		
# 		if sharedRecept == None:
# 			if RunName == 'SAMPL4':
# 				pathName,basename = os.path.split(path)
# 				pbits2 = pathName.split('/')
# 				receptName = pbits2[-2]
# #				 print receptName
# #				 receptorPath = pathName + '/' + receptName + '.pdbqt'
# 			else:
# 				print 'Local_Top_visit_ADV: missing receptor?!',RunName
# 				receptName = '??'
# 				continue
# 		else:
# 			rbits = sharedRecept.split('/')		  
# 			receptFile = rbits[-1]
# 			ppos = receptFile.find('.')
# 			receptName = receptFile[:ppos]

		dk = (exptname,receptName,ligand)

		allPaths[dk] = path
		relPath = path.replace(ProcDir,'')
		outs.write('%s,%s,%s,%s,%s\n' % (exptname,lib,receptName,ligand,relPath))
				
		if dk in dataTbl:
			print 'Local_Top_visit_ADV: dup dataKey?!',dk
			ndup += 1
			continue

		dataTbl[dk] = ligData
							   
		totParse += 1
	
	outs.close()		
	print 'Local_Top_visit_ADV: TotParse=%d NDup=%d' % (totParse,ndup)
	
	summf  = outdir+ ('ADV_summ_%07d.csv' % (1))
	interf = outdir+ ('ADV_inter_%07d.json' % (1))
	rptLocalData_ADV(dataTbl,summf,interf)
	

# arg string ala:
# 140905
# /mgl/storage/wcg/processed/ /mgl/storage/wcg/crawl/test/  --exptList "['Exp72']" --verbose
# 140906
# /mgl/storage/wcg/processed/ /mgl/storage/wcg/crawl/test/  --exptList "# arg string ala:
# 140905
# /export/wcg/processed/ /export/wcg/crawl/test/  --exptList "['Exp81','Exp82','Exp84']" --verbose

# 141030: local rikPR253
# /Data/sharedData/coevol-HIV/WCG/processed/rikPR253/processed /Data/sharedData/coevol-HIV/WCG/crawl/rikPR253 --verbose

# 160519
# --exptList "['Exp96-LEDGF']" --verbose

parser = argparse.ArgumentParser(description='crawl_ADV arguments')
# parser.add_argument('ADV_topDir',type=str,help='Path to crawling rootdir')
# parser.add_argument('outDir',type=str,help='Path to directory to contain result files')
parser.add_argument('--exptListStr', action="store",help='list of subset exptDir to crawl(string)')
parser.add_argument("--verbose",  action="store_true",help="increase output verbosity")
parser.add_argument("--recon",  action="store_true",help="Reconnaissance sniffing only")

### top-level run commands
if __name__ == '__main__':

	## mgl crawl of processed
	import socket
	HostName = socket.gethostname()
	if HostName == 'mgl0':
		print 'running on mgl0, good!'
		BaseDir = '/export/wcg/'
	
	elif HostName == 'mgl3':
		print 'running on mgl3, slow(:'
		BaseDir = '/mgl/storage/wcg/'
	
	elif HostName.startswith('hancock'):
		print 'running local on hancock'
		BaseDir = '/Data/sharedData/coevol-HIV/WCG/'

	elif HostName.startswith('mjq'):
		print 'running local on mjq'
		BaseDir = '/home/Data/coevol-HIV/WCG/'

	else:
		print 
		sys.exit( ('unknown host %s' % (HostName)) )

	args, unknown = parser.parse_known_args()
	if args.verbose:
		print 'crawl_ADV: arguments'
		# NB: args is a Namespace object; 
		argsDict = vars(args)
		for k,v in argsDict.items():
		 	print '\t%s = %s' % (k,v)
	
	if len(unknown)>0:
		print 'crawl_ADV: huh?! Unkown arguments=', unknown
		assert False # can't do break or return here!
	
	if args.exptListStr:
		exptList = eval(args.exptListStr)
	else:
		exptList = None
		 
	# print '# PROFILING run!!'
	# import cProfile
	# cProfile.run(('mglTop_visit_ADV(%s,verbose=%s)' % (exptList, args.verbose)), CrawlDir +'ADV_mgl0_profile.txt')

	ProcDir =  BaseDir + 'process2/'
	CrawlDir = BaseDir + 'crawl2/'
	config.VDWExcludedLigAtoms = ['H','C']
	dockfsuffix = '*_VS.pdbqt'
	tgzProc=False
	config.procVersion == '0.1'
	
	mglTop_visit_ADV(tgzProc, exptList, args.recon, args.verbose)
