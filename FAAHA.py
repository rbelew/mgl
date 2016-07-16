''' FAAHA: top-level full-workflow module for FAAH only
Created on 26 May 16

@author: rik
'''

 
from collections import defaultdict
import ConfigParser
import cPickle
import csv
import datetime
import glob
import json
import math	 
import os
import re
import socket
import string
import subprocess
import sys

import pwd
import grp
import os

import random
import numpy as np

import scipy.cluster.hierarchy as sch
# NB, using scipy's spatial distance (not SKLearn's metrics)
import scipy.spatial.distance as distance
from scipy.sparse.dok import dok_matrix
 
import matplotlib
from operator import irepeat
from sets import ImmutableSet
matplotlib.use('Agg') # ASSUME no windowing; prevent figures from popping up with savefig()
# Although many examples use pylab, it is no longer recommended.														 
# import pylab as p
import matplotlib.pyplot as p
 
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
 
import pandas as pd

import pybel
ob = pybel.ob
import recap3 as recap
import macrocycleopener
import fastcluster  


## other FAAHA modules
import config
import crawl_ADV
import plot
import getPDBQT

## Global variables, constants

config.RunName = ''
LowEDir = ''
InterTblDir = ''
RLIFDir = ''
L2FDir = ''
FragSimDir = ''
R2FDir = ''
LigCoordDir = ''
PlotDir = ''
ArffDir = ''
HIFDir = ''

InterTypes = ('hba', 'hbd', 'mtl','ppi','tpi','vdw')
BinaryITypes = ('hba', 'hbd', 'mtl','ppi','tpi')

FE_coeff_tors = 0.2983 # cf AD4.1_bound.dat

TooBigE = 30.
# 2DO: appropriate limits on LigEff, LigEnth?
TooBigLigEff = 30.
TooBigLigEnth = 30.

LibList = ['NF','EN','CB','AS','VM']

LigDOFTbl = {}


AADict = {'ASP':'D', 'GLU':'E', 'LYS':'K', 'HIS':'H', 'ARG':'R',
			   'GLN':'Q', 'ASN':'N', 'SER':'S', 'ASX':'B', 'GLX':'Z',
			   'PHE':'F', 'TRP':'W', 'TYR':'Y',
			   'GLY':'G', 'ALA':'A', 'ILE':'I', 'LEU':'L', 'CYS':'C',
			   'MET':'M', 'THR':'T', 'VAL':'V', 'PRO':'P',
			   'HID':'H', 'HIE':'H', 'HIP':'H',
			   'ASH':'D', 'GLH':'E',
			   'LYN':'K', 'ARN':'R',
			   'HOH':'U', 'CL': 'J' }

SmartsPat = None # initialized to ob.OBSmartsPattern() in bldLig2Frag()

## utilities
def basicStats(l):
	"Returns avg and stdev"
	if len(l) == 0:
		return(0.,0.)

	tot = 0
	for n in l:
		tot += n
	avg = float(tot) / len(l)

	sumDiffSq = 0.
	for n in l:
		sumDiffSq += (n-avg)*(n-avg)

	stdev = math.sqrt(sumDiffSq) / float(len(l))
	return (avg,stdev)

def jaccardSim(set1,set2):
	return float(len(set1.intersection(set2))) / len(set1.union(set2))

def freqHist(tbl):
	"Assuming values are frequencies, returns sorted list of (val,freq) items in descending freq order"
	def cmpd1(a,b):
		"decreasing order of frequencies"
		return cmp(b[1], a[1])

	
	flist = tbl.items()
	flist.sort(cmpd1)
	return flist

def entropy2(ix,iy):
	x = float(ix)
	y = float(iy)
	if x==0:
		e=0.
	else:
		e =  -(x / (x+y)) * math.log( x / (x+y) , 2.0)
	if y > 0.:
		e += -(y / (x+y)) * math.log( y / (x+y) , 2.0)
	return e

def infoGain(tp,fp,tn,fn):
	all = tp+fp+tn+fn
	prob = float(tp+fp) / all
	pnot = 1. - prob
	ig = entropy2(tp+fn,fp+tn)
	pos = entropy2(tp,fp)
	neg = entropy2(fn,tn)
	ig -= prob * pos
	ig -= pnot * neg
	
	assert ig >= 0., 'Negative info gain?!'
	
	return (ig,prob,pos,neg)

def ranges2list(rangeList):
	'invert bldRanges()'
	l = []
	for e in rangeList:
		if type(e)==type(1): # singleton int
			l.append(e)
		else:
			l += range(e[0],e[1]+1)
	return l

def rndEstr(e,ndig=0):
	"round to nearest integral +ndig kcal"
	fstr = '{0:0%d.%df}' % (ndig+4,ndig) # so two-digit energies sort correctly as strings!
	rstr = fstr.format(e)
	return rstr

def touchTimeDiff(f1,f2):
	'''given two files (earlyFile, lateFile), return seconds between their touch times
	'''
	t1 = os.stat(f1).st_mtime
	t2 = os.stat(f2).st_mtime
	
	return t2 - t1

def setConfigOptions(cfg,module='',verbose=True):
	'''bind all config options to variables (in module)
	option type (Pred,Int,Float) stripped from variable name before binding
	ASSUME cfg has cfg.optionxform=str so case is maintained
	'''

	# NB: DEFAULT parameters included by ALL sections!
	if verbose: print '<setConfigOptions>'
	for section in cfg.sections():
		for k,v in cfg.items(section):
			# import pdb; pdb.set_trace()
			if v == 'None':
				v = None
			elif k.endswith('Pred'):
				v = cfg.getboolean(section,k)
				k = k[:-4]
			elif k.endswith('Int'):
				v = cfg.getint(section,k)
				k = k[:-3]
			elif k.endswith('Float'):
				v = cfg.getfloat(section,k)
				k = k[:-5]
			else:
				v = '"' + v + '"'
			# print '%s %s' % (k,v),
			if module=='':
				exec('%s = %s' % (k,v))
				if verbose: print '%s \t = %s' % (k,eval(k))
			else:
				exec('%s.%s = %s' % (module,k,v))
				if verbose: print '%s.%s \t = %s' % (module,k,eval('%s.%s' % (module,k)))
				
	if verbose: print '</setConfigOptions>'

		
## eo utilities	
	
def bldExptTbl(inf):
	'''return summary table of experiments
	
	exptKey = (exptNo,prot,recept,site,lib)
	experiment attributes: bstart bend sys lib protein receptor site
	'''
	
	reader = csv.DictReader(open(inf))
	# bstart,bend,receptor,protein,site,library,experiment,pdb,gpf,dpf,program,config
	# 5398,5398,x1k6_EqMD,PR,Exo_ChnA,CB_bb,27,1KZK,faah_ExpandedExo_4_Jan2009.gpf,faah_template_EndByGen_March09.dpf,AD

	redunKeyFnd = False
	exptTbl = {}
	for i,entry in enumerate(reader):
		
		if entry['bstart'].startswith('# '):
			continue
		
		# HACK:	  
		# replace any '_' in receptor, site, library with '-'
		# to avoid interactions with underbar concatenations later
		
		exptData = {}

		try:
			exptNo = entry['experiment']
			# 160217: require consistent hyphens within prot (experiment) name
			# prot = entry['protein'].replace('_','-')
			prot = entry['protein']
			if prot.find('_') != -1:
				sys.exit('bldExptTbl: require hyphens within protein name?!')
				
			recept = entry['receptor'].replace('_','-')
			recept = entry['receptor'].replace('.pdbqt','')
			site = entry['site'].replace('_','-')
			lib = entry['library'].replace('_','-')
			
			exptData['bstart']  = int(entry['bstart'])
			exptData['bend']	= int(entry['bend'])
			exptData['sys']	 = config.RunType
		except Exception,e:
			print 'bldExptTbl: bad entry?! %s %s' % (e,entry)
			continue
		
		exptData['lib']	 = lib
		exptData['protein'] = prot
		exptData['receptor'] = recept
		exptData['site'] = site

		exptKey = (exptNo,prot,recept,site,lib)
		if exptKey in exptTbl:
			print 'bldExptTbl: dup key?! %s \n\twas: %s\n\tnew: %s' % (exptKey,exptTbl[exptKey],exptData)
			redunKeyFnd = True
			continue
		
		exptTbl[exptKey] = exptData
		
	print 'bldExptTbl: NExpt=%d' % (len(exptTbl))
	if redunKeyFnd:
		print 'bldExptTbl: redun keys disallowed!'
		return None
	
	return exptTbl
	
def FAAHA_expt(exptTbl,exptList,faahDirParent,outDir):
	''' Identify top ncand and top frac4Thresh for each experiment
	
		first pass sorts all energies, identifies threshold, 
		outputs best ncand  candidates ACROSS ALL EXPERIMENTS to best_*.csv
		two thresholds computed: threshold=thresh@frac4Thresh; bestCandThresh=thresh@ncand
		computes, saves NonZincLigTbl in getThresh()
		NB: exptTbl augmented with threshold!
	'''
				
	thresh = config.Frac4Thresh
	ncand = config.NCand
	dcrit = config.Criterion
		
	## Process experiments in exptList
		
	# print 'Expt,Lib,NBatch,NRcd,NLigand,NLowE,NDup,AvgRcd,FracLowE' 
			
	allExpt = exptTbl.keys()
	expt2do = []
	for exptKey in allExpt:
		(exptNo,prot,recept,site,lib) = exptKey
		# NB: accept exptList with 'Exp' ala processed, crawled directories
		if exptList:
			if exptNo in exptList:
				expt2do.append(exptKey)
			else:
				print 'skipping',exptKey
		else:
			expt2do.append(exptKey)
			
	print 'FAAHA_expt: exptList=%s NExpt=%d' % (exptList,len(expt2do))

	candf = outDir+('bestLig_%d.csv' % ncand)

	cands = open(candf,'w')
	cands.write('Expt,Batch,Ligand,E,FullExpt\n')
	
	runType = config.RunType
	newExptTbl = {}
	expt2do.sort()
	for exptKey in expt2do:
		exptNo,prot,recept,site,lib = exptKey
		exptData = exptTbl[exptKey]
		exptName = bldExptStr(exptKey)

		# Experiments collapse all results into single file; 
		batchList = ranges2list( [(exptData['bstart'],exptData['bend'])] )			
			
		lib = exptTbl[exptKey]['lib']

		nbatch = len(batchList)
		
		faahDir = faahDirParent + exptNo + '/' + recept + '/'
					
		frac4Thresh=config.Frac4Thresh
		uniqLigTbl, thresh, bestCandThresh = getThresh(exptKey,runType,batchList,faahDir,frac4Thresh, \
													   ncand=config.NCand,dcrit=config.Criterion,initRun=True)
				
		newExptTbl[exptKey] = exptData.copy()
		
		# NB: exptTbl augmented with threshold!
		newExptTbl[exptKey]['thresh'] = thresh
		# 150930  keep ncand thresh info too
		newExptTbl[exptKey]['ncand'] = ncand
		newExptTbl[exptKey]['bestCandThresh'] = bestCandThresh

		## 2d pass

		outs = open(LowEDir+('%s_lowE.csv' % exptName),'w')
		# cf. mglTop_visit_AD() for header line
		fldNameLine = 'Batch,Ligand,E,Eff,DCVal,Nvdw,Ninter\n'
		outs.write(fldNameLine)
	
		nrcd = 0
		nlowE = 0
		ndup = 0
		nbadE = 0
		nmissUniq = 0
		nbadZinc = 0
		
		lowETbl = {}
		ligTblE = {}
		runType = config.RunType
				
		exptLigDir = EliteDir + exptName + '/'

		nnrtiDrop = 0
			
		for isf,batchNo in enumerate(batchList):
			summPath = faahDir+('summ/%s_summ_%07d.csv' % (runType,batchNo))
				
			try:
				inStr = open(summPath)
			except Exception,e:
				# ASSUME counted as nmissf above
				# print "FAAHA_expt: can't open2 %s?!" % (summPath)
				continue
		 
			for il,line in enumerate(inStr.readlines()):
				if il == 0:
					# Expt,Recept,Ligand,E,Eff,Nvdw,Ninter
					# written by get_ADInfo.rptData()
					continue
				nrcd += 1
				try:
					flds =line[:-1].split(',')
					(expt,batch,ligRecept,ligand,e,eff,nvdv,ninter) = flds

					if dcrit == 'energy':
						dval = float(e)
					elif dcrit == 'ligEff':
						dval = float(eff)
					elif dcrit == 'ligEnth':
						# Assume LigDOFTbl already loaded
						if ligand not in LigDOFTbl:
							print 'FAAHA_expt: lig missing from LigDOFTbl?!',ligand
							continue
						dof = LigDOFTbl[ligand]
						natom = float(e) / float(eff)
						natomi = int(natom)
						assert natomi-natom > 1e-2, 'huh?'
							
						dval = (float(e) + (dof * FE_coeff_tors)) / natom

				except Exception, exp:
					print "FAAHA_expt: bad line(2) %s %d %s?!" % (summPath,il,exp)
					break # don't try to read other lines from inStr
				
				# Pass2: need to make this thresh criterion-sensitive, too
				
				if (dcrit == 'energy' and abs(dval) > TooBigE) or \
					(dcrit == 'ligEff' and abs(dval) > TooBigLigEff) or \
					(dcrit == 'ligEnth' and abs(dval) > TooBigLigEnth):
					nbadE += 1
					continue
								
				ligIdx = normLigand(ligand)
				if ligIdx == -1:
					nbadZinc += 1
					continue
							
				ligand2 = ligIdx2zinc(ligIdx)
				
				if ligand2 != ligand:
					# NB: two ligand cleanups!
					if not (ligand.endswith('.VS') or ligand.find('pras') != -1):
						print 'FAAHA_expt: bad ligand indexing?!', ligand, ligIdx, ligand2
						import pdb; pdb.set_trace()
						continue

				# 160122: Remove redundant ligands
				if ligIdx not in uniqLigTbl:
					nmissUniq += 1
					print 'FAAHA_expt: ligand not seen by getThresh()?!',ligand, batchNo
					continue
				
				# only keep one (first) docking result wrt/ same ligand
				ubatch,udval = uniqLigTbl[ligIdx]
				if ubatch != batchNo:
					ndup += 1
					continue
	
				## thresh
				# NB: thresh set EQUAL to highest energy in getThresh()
				if dval <= thresh:
					nlowE += 1
					## NB: all 3 vals (e,eff,dval) written to lowE file; why not!
					newline = '%d,%s,%s,%s,%s,%s,%s\n' % (batchNo,ligand,e,eff,dval,nvdv,ninter)	  
					outs.write(newline)
				if dval < bestCandThresh:
					cands.write('%s,%d,%s,%f,%s\n' % (exptNo,batchNo,ligand,dval,exptName))
				
			inStr.close() # eo-batch file
			
		nlig = len(uniqLigTbl)
			
		outs.close()
					   
		try:
			avgRcd = float(nrcd)/nbatch
		except:
			avgRcd = 0
						
		print 'FAAHA_expt: %s NBatch=%d NRcd=%d NLigand=%d NBadZinc=%d NNonUniq=%d NLowE=%d NDup=%d NBadE=%d AvgRcd/File=%f NNRTIDrop=%d'  % \
				(exptName,nbatch,nrcd,nlig,nbadZinc,nmissUniq,nlowE,ndup,nbadE,avgRcd,nnrtiDrop)
		
	cands.close()
	return newExptTbl

def getThresh(exptKey,runType,batchList,faahDir,frac4Thresh,ncand,initRun=False,dcrit='energy'):
	'''set threshold, return uniqLigTbl, 
	141219: selection criterion made a parameter
	'''
	
	exptNo,prot,exptRecept,site,lib = exptKey

	exptName = bldExptStr(exptKey)

	allE = []
	noddE = 0
	nmissf = 0
	nrcd=0
	ndup=0
	ndupOut=0
	nnrtiDrop = 0		
	nbadZinc = 0
	
	ligTblE = {} # ligidx -> [batchNo,dval]

	if initRun:
		dupLigFile = SummRptDir + 'dupDockedLig.csv'
		# multiple experiments share same dupLigFile
		if os.path.isfile(dupLigFile):
			dupStream = open(dupLigFile,'a')
		else:
			dupStream = open(dupLigFile,'w')
			dupStream.write('Ligand, PrevBNo, PrevDval, BNo, Dval\n')

		config.NonZincLigTbl = {} # ligand name -> nonzincID 
		config.NZIdx2LigTbl = {} # nonzincID -> ligand name
		config.NNonZincLig = 0
		
	else:
		loadNonZinc(exptName)
	
	for isf,batchNo in enumerate(batchList):
		summPath = faahDir+('summ/%s_summ_%07d.csv' % (runType,batchNo))
			
		try:
			inStr = open(summPath)
		except Exception,e:
			print "getThresh: can't open1 %s?!" % (summPath)
			nmissf += 1
			continue
		
		for il,line in enumerate(inStr.readlines()):
			if il == 0:
				# Expt,Recept,Ligand,E,Eff,Nvdw,Ninter
				# written by crawl_ADV.rptData_ADV()
				continue

			nrcd += 1
			try:
				flds =line[:-1].split(',')
				(expt,batch,ligRecept,ligand,e,eff,nvdv,ninter) = flds
				
				## 141219: Pass 1: selection criterion made a parameter
				
				if dcrit == 'energy':
					dval = float(e)
				elif dcrit == 'ligEff':
					dval = float(eff)
				elif dcrit == 'ligEnth':
					# Assume LigDOFTbl already loaded
					if ligand not in LigDOFTbl:
						print 'getThresh: lig missing from LigDOFTbl?!',ligand
						continue
					dof = LigDOFTbl[ligand]
					natom = float(e) / float(eff)
					natomi = int(natom)
					assert natomi-natom > 1e-2, 'huh?'
						
					dval = (float(e) + (dof * FE_coeff_tors)) / natom	
											
			except Exception, e:
				print "getThresh: bad line(1) %s %d %s?!" % (summPath,il,e)
				noddE += 1
				break # don't try to read other lines from inStr
			
			# 160122: Remove redundant ligands
			ligIdx = normLigand(ligand)
			if ligIdx == -1:
				nbadZinc += 1
				continue
			
			if ligIdx in ligTblE:
				ndup += 1
				if ligTblE[ligIdx][1] != dval and initRun:
					# print 'getThresh: same ligand with differing dval?!',ligand, ligTblE[ligIdx][0], ligTblE[ligIdx][1], batchNo, dval
					dupStream.write('%s,%s,%s,%s,%f\n' % (ligand, ligTblE[ligIdx][0], ligTblE[ligIdx][1], batchNo, dval))
					ndupOut += 1
				continue
			else:
				ligTblE[ligIdx] = [batchNo,dval]

			# Pass1: need to make this thresh criterion-sensitive, too
			
			if (dcrit == 'energy' and abs(dval) > TooBigE) or \
				(dcrit == 'ligEff' and abs(dval) > TooBigLigEff) or \
				(dcrit == 'ligEnth' and abs(dval) > TooBigLigEnth):
				noddE += 1
				continue
								
			allE.append(dval)
			
		inStr.close() # eo-batch file
	
	if initRun:
		dupStream.close() # eo all batches

		if config.NNonZincLig>0:
			# cf. bldNonZincIdx()
				
			allLig = config.NonZincLigTbl.keys()
			allLig.sort()
			
			nzligFile = LowEDir + 'nonZincLig.csv'
			print 'getThresh: Writing %d NonZinc ligands to %s' % (config.NNonZincLig,nzligFile)
			
			# multiple experiments share same nzligFile
			if os.path.isfile(nzligFile):
				outs = open(nzligFile,'a')
			else:
				outs = open(dupLigFile,'w')
				outs.write('Ligand,NZIdx\n')
							
			for lig in allLig:
				outs.write('%s,%d\n' % (lig,config.NonZincLigTbl[lig]))
			outs.close()

	# NB: dval used as generic value from here forward
	# "E" name is no longer correct, but doesn't hurt anything!

	if len(allE)==0:
		# print 'getThresh: No good summ files1?! Expt=%s NBatch=%d' % (exptName,len(batchList))
		# return ({}, 0.,0.)
		sys.exit('getThresh: No good summ files1?! Expt=%s NBatch=%d' % (exptName,len(batchList)))
	
	# ASSUME all algorithms require n log n??
	allE.sort()
	if frac4Thresh==1.0:
		# NB: thresh set to EQUAL highest energy found
		thresh = allE[-1]
		threshIdx = len(allE)
	else:
		threshIdx = int(round(float(len(allE) * frac4Thresh)))
		thresh = allE[threshIdx]	   
	
	if ncand > len(allE):
		print 'getThresh: ncand < allLig=%d; using all' % (len(allE))
		ncand = len(allE)-1
	bestCandThresh = allE[ncand]

	print 'getThresh: Expt=%s NLigand=%d NBadZinc=%d NDup=%d NDupOut=%d NMissf=%d NOddE=%d ThreshFrac=%5.2f ThreshIdx=%d Thresh=%f BestCandThresh=%f NNRTIDrop=%d' % \
					  (exptName,len(ligTblE),nbadZinc,ndup,ndupOut,nmissf,noddE,frac4Thresh,threshIdx,thresh,bestCandThresh,nnrtiDrop)
	
	return ligTblE, thresh, bestCandThresh
	
def bldZincIdx(zstr):
	# ASSUME no more than two digit ZINC suffices
	if zstr.find('_') != -1:
		bpos = zstr.find('_')
		zno = int(zstr[4:bpos])
		# 2do where is "ZINC06648739_10_" with trailing '_' getting produced?!
		zsuf = zstr[bpos+1:].strip(' _')  
		if len(zsuf)==0:
			zsuf = 0
		else:
			zsuf = int(zsuf)
		zidx = 100 * zno + zsuf
	else:
		zno = int(zstr[4:].strip())
		zidx = 100 * zno
	return zidx

def ligIdx2zinc(ligIdx):
	zid = ligIdx / 100
	sufId = ligIdx % 100

	if ligIdx < 0:
		if ligIdx in config.NZIdx2LigTbl:
			return config.NZIdx2LigTbl[ligIdx]
		else:
			print 'ligIdx2zinc: missing nonZincIdx?!',ligIdx
			return 'Lig??'

	zincid = 'ZINC%08d' % (zid)

	if sufId != 0:
		zincid += '_%d' % (sufId)

	return zincid

def bldNonZincIdx(lig):
	'assign arbitrary (negative) numbers to non-Zinc ligands'
	if lig in config.NonZincLigTbl:
		return config.NonZincLigTbl[lig]
	else:
		config.NNonZincLig += 1
		config.NonZincLigTbl[lig] = (-(config.NNonZincLig))
		config.NZIdx2LigTbl[(-(config.NNonZincLig))] = lig
		return  (-(config.NNonZincLig))
	
def normLigand(lig):
	'''HACK: reduce ZINC id to integer, for reduced storage
	'''
	
	if lig.find('ZINC') != -1:
		try:
			zincIdx = bldZincIdx(lig)
		except:
			print 'normLigand: odd ligand?!',lig
			# import pdb; pdb.set_trace()
			return -1
	else:
		# 160526: no nonZINC allowed!
		# zincIdx = bldNonZincIdx(lig)
		return -1
		
	return zincIdx
	
def getIsoRAtom(raa,ratom):
	'''map any references to (ratom2) -> (ratom1) for 
	OD2 in D, OE2 in E, and CD2 or CE2 in F and Y
	'''
	
	if raa == 'D':
		if raa == 'OD2':
			isoRAtom = 'OD1'
	elif raa == 'E':
		if raa == 'OE2':
			isoRAtom = 'OE1'
	elif raa == 'F' or raa == 'Y':
		if raa == 'CD2':
			isoRAtom = 'CD1'
		elif raa == 'CE2':
			isoRAtom = 'CE1'
	else:
		isoRAtom = ratom
		
	return isoRAtom
	
	
def bldFocalInterTbl(exptName,faahDir,batch2ligTbl):
	'''retrieval of interaction data for just focal ligands (collected by batch)
	returns receptInterTbl: (rchain,raa,ratom) -> {itype - > {liname -> [ ligID ] } }
	NB: receptor AA naming normalized 
	NB: ALL interactions included in interTbl; 
		exclusions applied later via bldFeature2() in analRLIF(), analBldHIFeature
	'''
	
	# NB: receptInterTbl retains lig atom distinctions
	receptInterTbl = defaultdict( lambda: defaultdict(lambda: defaultdict( list ))) 
	# (rchain,raa,ratom) -> {itype - > {liname -> [ ligID ] } }
		 
	nrcd = 0
	nmissf = 0
	nactive = 0
	nVDWonly = 0
	ligFndSet = set()
	
	for ib,bno in enumerate(batch2ligTbl.keys()):
		
		inf = faahDir+('inter/%s_inter_%07d.json' % (config.RunType,int(bno)))
			
		try:
			inStr = open(inf)
			allInter = json.load(inStr)
			inStr.close()
		except:
			nmissf += 1
			continue
		
		focalLigTbl = batch2ligTbl[bno]
		
		if len(focalLigTbl)==0:
			print 'bldFocalInterTbl: empty batch?!',ib,bno
			continue
		
		for il,interInfo in enumerate(allInter):
			nrcd += 1
			
			# cf crawl_ADV.rptData_ADV()
			# [ [Expt,BatchNo,Recept,Lig, [IType,[InterEnum] ] ] ]
			(expt,batchNo,recept,ligand, interList) = interInfo
			
			ligIdx = normLigand(ligand)   
								 
			if ligIdx not in focalLigTbl:
				continue

			nonVDWfnd = False
			for interTypeList in interList:
				itypeIdx, interEnum = interTypeList
				itype = InterTypes[itypeIdx]
				if itype != 'vdw':
					nonVDWfnd = True
				for iinfo in interEnum:
					if itype == 'vdw':
						if config.vdwLigandAtom:
							(rchain,raa,ratom,liname) = iinfo
						elif config.ADFeatures == 'nonVDW':
							continue
						else:
							(rchain,raa,ratom) = iinfo
							liname = ''
					elif itype == 'tpi' or itype=='ppi':
						(rchain,raa,ratom,liname) = iinfo
						# NB: drop ligand ring center; use X's to fill out feature
						ratom = 'X'
						liname = 'X'
					else:
						(rchain,raa,ratom,liname) = iinfo

					# NB: normalize receptor AA naming
					raas = bldFeaturePrefix(rchain, raa)
					rchain, newRAA = raas.split('_')
					
					## 160610: Handle isomorphic atom naming
					# References made to EITHER of alternative ratom names treated as reference to FIRST NAME
					# NB: no checking here for interactions with BOTH of ambiguous ratom names
					
					if config.MapIsoRAtoms:
						isoRAtom = getIsoRAtom(raa,newRAA)
					else:
						isoRAtom = newRAA
					 
					k = (rchain,newRAA,isoRAtom)
		
					receptInterTbl[k][itype][liname].append( ligIdx )
					ligFndSet.add(ligIdx)
						
			if not nonVDWfnd:
				nVDWonly += 1

	print 'bldFocalInterTbl: %s NLigand=%d NVDWOnly=%d N(Chain+RAA+RAtom)=%d NMissFile=%d NAllLig=%d' % \
		(exptName,len(ligFndSet),nVDWonly,len(receptInterTbl),nmissf,nrcd)

	## NB: need to make serializable, for pickel!
	
	interTbl2 = {}
	for rfeat,ligTbl in receptInterTbl.items():		 
		# (rchain,raa,ratom) -> {itype - > {latom -> [ ligID ] } }
		newLigTbl = {}
		for itype,ligFeatTbl in ligTbl.items():
			newLFTbl = {}
			for latom,ligIDList in ligFeatTbl.items():
				newLFTbl[latom] = ligIDList[:]
#				 for ligIdx in ligIDList:
#					 if ligIdx<0 and itype != 'vdw':
#						 print '2,%s,%s,%s,%s' % (rfeat,itype,latom,ligIdx)
					
				
			newLigTbl[itype] = newLFTbl
			assert len(newLFTbl)==len(ligFeatTbl), 'bldFocalInterTbl: bad newLFTbl?!'
		assert len(newLigTbl)==len(ligTbl), 'bldFocalInterTbl: bad newLigTbl?!'
		interTbl2[rfeat] = newLigTbl
	assert len(interTbl2)==len(receptInterTbl), 'bldFocalInterTbl: bad interTbl2?!'
	
	return interTbl2

def analBldHIFeatures(faahDir,exptName,bnoList,thresh,r2lTbl,uniqLigTbl,featFile):	
	'''build energy discrimination predicate
	v3 use thresh, dont build explicit energy distribution
	160122: summ files may contain redundant references to ligands!
	use classTbl to maintain ALL referenced ligands (vs nrcd)
	'''
	
	print "analBldHIFeatures: %s NBatch=%d" % (exptName,len(bnoList))

	classTbl = {} # ligand -> aboveThreshP
	noddE = 0
	npos = 0
	nmissf = 0
	ndup = 0
	ncdup = 0
	nmissUniq = 0
		
	nnrtiDrop = 0

	## first pass: build energy distribution
	
	for ib,bno in enumerate(bnoList):

		summPath = faahDir+('summ/%s_summ_%07d.csv' % (config.RunType,bno))

		try:
			inStr = open(summPath)
		except:
			nmissf += 1
			# import pdb; pdb.set_trace()
			continue
		
		for il,line in enumerate(inStr.readlines()):
			if il == 0:
				continue
			flds =line[:-1].split(',')

			try:
				(expt,batch,recept,ligand,e,eff,nvdv,ninter) = flds
			except Exception,e:
				print "analBldHIFeatures: bad line %s %d %s?!" % (summPath,il,e)
				noddE += 1
				break # don't try to read other lines from inStr
				
			# HACK: reduce ZINC id to integer, for reduced storage
			
			zincIdx = normLigand(ligand)
						
			if abs(float(e)) > TooBigE:
				# print 'FAAHA: odd energy?!',isf,summF,il,e
				noddE += 1
				continue

			# 160122: Remove redundant ligands
			if zincIdx not in uniqLigTbl:
				nmissUniq += 1
				# print 'analBldHIFeatures: ligand not seen by getThresh()?!',ligand, bno
				continue
			
			# only keep one (first) docking result wrt/ same ligand
			ubatch,udval = uniqLigTbl[zincIdx]
			if ubatch != bno:
				ndup += 1
				continue
					
			if zincIdx in classTbl:
				ncdup += 1
				tst = float(e)<=thresh
				if tst != classTbl[zincIdx]:
					print "analBldHIFeatures: differing test for same ligand?!",ligand,classTbl[zincIdx]
				continue
				
			classTbl[zincIdx] = float(e)<=thresh
			if float(e)<=thresh:
				npos += 1
					   
		inStr.close()

	print "analBldHIFeatures: %s NLigand=%d NMissLig=%d NDupLig=%d NDupCat=%d NOddE=%d Thresh = %f NPos=%d NMissFile=%d NDropNRTI=%d" % \
		(exptName,len(classTbl),nmissUniq,ndup,ncdup,noddE,thresh,npos,nmissf,nnrtiDrop)	
		
	if len(classTbl)==0:
		print 'analBldHIFeatures: no ligands?!'
		return None

	featInfoTbl = {}
	ncnt = 0
		
	for rlif,ligList in r2lTbl.items():			
		#  rlif -> [(ligIdx,latom)]
		if rlif not in featInfoTbl:
			featInfoTbl[rlif] = [0,0]
				 
		for ligInfo in ligList:
			ligIdx,latom = ligInfo

			ncnt += 1
				
			if classTbl[ligIdx]:
				featInfoTbl[rlif][0] += 1
			else:
				featInfoTbl[rlif][1] += 1

	print "analBldHIFeatures: %s NLigand=%d NPos=%d NRLIFxLig=%d Nfeatures=%d" % \
		(exptName,len(classTbl),npos,ncnt,len(featInfoTbl))
		
	allFeat = featInfoTbl.keys()
	allFeat.sort()
	
	nbadEntropy = 0
	totLig = len(classTbl)
	
	# Read by analCluster()
	outs = open(featFile,'w')
	outs.write('F,Pr,Info,PosEnt,NegEnt\n')
	for f in allFeat:
		# TP: in low-mode,	 feature present
		# FN: in low-mode,	 no feature
		# FP: not in low-mode, feature present
		# TN: not in low-mode, no feature
		tp = featInfoTbl[f][0]
		fn = npos - tp
		fp = featInfoTbl[f][1]
		tn = totLig - npos - fp
		
		all = tp+fp+tn+fn
#		 assert (fn>0 and tn>0), 'analBldHIFeatures: bad counts1?!'
#		 assert all == nrcd, 'analBldHIFeatures: bad counts2?!'
		if not(fn>=0 and tn>=0 and all == totLig):
			nbadEntropy += 1
			outs.write('"%s",,,,,?1,%d,%d,%d,%d\n' % (f,tp,fp,npos,totLig))
			continue

		try:
			ig, prFeat, posEnt,negEnt = infoGain(tp,fp,tn,fn)
		except Exception,ex:
			# print 'analHIFeatures: badEntropy?!',exptName,ex,ifeat,npos,nrcd,tp,fp,tn,fn
			nbadEntropy += 1
			outs.write('"%s",,,,,?2,%d,%d,%d,%d\n' % (f,tp,fp,npos,totLig))
			continue
					
		outs.write('"%s",%f,%f,%f,%f\n' % \
				   (f,prFeat,ig,posEnt,negEnt))
	outs.close()
	print "analBldHIFeatures: %s NBadEntropy=%d" % (exptName,nbadEntropy)

def loadRLIFInfo(hifFile):
	'''return rlifInfoTbl: rlif -> info
	'''
	nerr = 0
	rlifTbl = {}
	reader = csv.DictReader(open(hifFile))
	# first pass: read them all
	for i,entry in enumerate(reader):
		# F,Pr,Info,PosEnt,NegEnt
		
		if None in entry:
			# cf analBldHIFeatures() error reporting
			print 'loadRLIFInfo: bad info?!',entry['F'],entry[None]
			nerr += 1
			continue
		
		rlifTbl[ entry['F'] ] = float( entry['Info'] )

	if nerr>0:
		print 'loadRLIFInfo: NErr=%d %s' % (nerr,hifFile)
	return rlifTbl

def analExptTblRLIFTypeInfo(exptTbl,HIFDir,outf):

	allExpt = exptTbl.keys()
	allExpt.sort()

	allRLIFTbl = {} # expt -> RLIF -> info
	for expt in allExpt:
		
		exptNo,prot,recept,site,lib = expt
		exptName = bldExptStr(expt)
	
		hifFile = HIFDir+('%s.csv' % (exptName))
		rlifTbl = loadRLIFInfo(hifFile)
		allRLIFTbl[exptName] = rlifTbl
	
	infoBuckets = []
	for dp in range(1,7):
		for m in [4.,2.,1.]:
			infoBuckets.append( (m * ((0.1) ** dp)) )
	infoBuckets.append(0.)
	
	rlifStats = {} # expt -> itype -> countVec
	for exptName in allRLIFTbl.keys():
		rlifStats[exptName] = {}
		for itype in InterTypes:
			rlifStats[exptName][itype] = [0 for v in infoBuckets]
			
	for exptName in allRLIFTbl.keys():
		for rlif in allRLIFTbl[exptName]:
			info = allRLIFTbl[exptName][rlif]
			bits = rlif.split('_')
			itype = bits[3]
			for i,v in enumerate(infoBuckets):
				if info > v:
					rlifStats[exptName][itype][i] += 1
					break
	
	outs = open(outf,'w')
	outs.write('Info')
	for exptName in allRLIFTbl.keys():
		for itype in InterTypes:
			outs.write(',%s_%s' % (exptName,itype))
	outs.write('\n')
	for i,v in enumerate(infoBuckets):
		outs.write('%e' % v)
		for exptName in allRLIFTbl.keys():
			for itype in InterTypes:
				outs.write(',%d' % rlifStats[exptName][itype][i])
		outs.write('\n')
	outs.close()		 

def analListRLIFITypeInfo(exptList,outf):

	allRLIFTbl = {} # expt -> RLIF -> info
	for exptName,hifFile in exptList:
			
		rlifTbl = loadRLIFInfo(hifFile)
		allRLIFTbl[exptName] = rlifTbl
	
	infoBuckets = []
	for dp in range(1,7):
		for m in [4.,2.,1.]:
			infoBuckets.append( (m * ((0.1) ** dp)) )
	infoBuckets.append(0.)
	
	rlifStats = {} # expt -> itype -> countVec
	for exptName in allRLIFTbl.keys():
		rlifStats[exptName] = {}
		for itype in InterTypes:
			rlifStats[exptName][itype] = [0 for v in infoBuckets]
			
	for exptName in allRLIFTbl.keys():
		for rlif in allRLIFTbl[exptName]:
			info = allRLIFTbl[exptName][rlif]
			bits = rlif.split('_')
			itype = bits[3]
			for i,v in enumerate(infoBuckets):
				if info > v:
					rlifStats[exptName][itype][i] += 1
					break
	
	outs = open(outf,'w')
	outs.write('Info')
	for exptName in allRLIFTbl.keys():
		for itype in InterTypes:
			outs.write(',%s_%s' % (exptName,itype))
	outs.write('\n')
	for i,v in enumerate(infoBuckets):
		outs.write('%e' % v)
		for exptName in allRLIFTbl.keys():
			for itype in InterTypes:
				outs.write(',%d' % rlifStats[exptName][itype][i])
		outs.write('\n')
	outs.close()		 
		
	  
def loadHInfoRLIF(hifFile,maxnhif=-1,minInfo=-1.):
	'''return all RLIF with info > minInfo or first maxnhif
	'''
	rlifTbl = {}
	reader = csv.DictReader(open(hifFile))
	# first pass: read them all
	for i,entry in enumerate(reader):
		# F,Pr,Info,PosEnt,NegEnt
		rlifTbl[ entry['F'] ] = float( entry['Info'] )
	
	nrlif = len(rlifTbl)
	allRLIF = rlifTbl.keys()
	allRLIF.sort(key=lambda k: rlifTbl[k], reverse=True)
	
	if maxnhif==-1 and minInfo==-1.:
		print 'loadHInfoRLIF: returning all %d' % (len(allRLIF))
		return allRLIF
	
	if maxnhif > -1:
		print 'loadHInfoRLIF: returning first %d/%d' % (maxnhif,len(allRLIF))
		return allRLIF[:maxnhif]
	
	minRLIFIdx = None
	for ir,rlif in enumerate(allRLIF):
		if rlifTbl[rlif] < minInfo:
			minRLIFIdx = ir
			break
	print 'loadHInfoRLIF: returning info>%f %d/%d' % (minInfo,minRLIFIdx,len(allRLIF))
	return allRLIF[:minRLIFIdx]
	
def bldExptStr(exptKey):
	s = '_'.join(exptKey)
	return s 

def analRLIF(exptKey):
	'''Create intertbl, RLIF2Ligand and RLIFFreq  
	'''

	# faahDirParent -> CrawlDir
	# outDir -> SummRptDir

	exptNo,prot,recept,site,lib = exptKey
	exptName = bldExptStr(exptKey)
	
	r2lf = InterTblDir + '%s_r2l.csv' % (exptName)
	
	# ASSUME existence of r2lf file implies all vdwPatch, etc also constructed
	if os.path.isfile(r2lf):
		print 'analRLIF: %s r2lf exists; skipping' % (exptName)
		return
	
	faahDir = CrawlDir +  '%s/%s/' % (exptNo,recept)

	ligTbl, batch2ligTbl = provideEliteDist(exptName)

	itpklf = InterTblDir + '%s_interTbl.pkl' % (exptName)
	if os.path.isfile(itpklf):
		print 'analRLIF: %s interTbl exists, loading' % (exptName)
		interTbl = cPickle.load(open(itpklf,'rb'))
		
	else:
		print 'analRLIF: %s interTbl does not exist, building' % (exptName)
		
#			interTbl = bldFocalInterTbl(faahDir,batch2ligTbl,activeLig)
		interTbl = bldFocalInterTbl(exptName,faahDir,batch2ligTbl)
		cPickle.dump(interTbl, open(itpklf,'wb'))

	# interTbl: (rchain,raa,ratom) -> {itype - > {liname -> [ ligID ] } }

	allRA = set()
	# Initialize rlifTbl here
	rlifTbl = defaultdict(list) # RLIF -> [ligIdx]

	if not config.vdwLigandAtom:
		print 'analRLIF: vdwLigandAtom=False; not forming vdwPatches'
	elif config.ActiveTargetLib:
		print 'analRLIF: ActiveTargetLib; not forming vdwPatches'
	else:
		la2vdwTbl = bldLigAtom2vdw(interTbl) # ligIdx -> latomFull -> [(rchain,raa,ratom)]
		  
		vdwpTbl = {} # vdwPatchTbl: raSet -> latype -> freq
		for il,ligIdx in enumerate(la2vdwTbl.keys()):
			la2raTbl = defaultdict(set) # ligAtom -> set(ratom)
			# first pass thru all ligatoms in ligIdx: collect all ratoms mentioning it
			for latomFull in la2vdwTbl[ligIdx]:
				for ratomFull in la2vdwTbl[ligIdx][latomFull]:
					rkeys = '_'.join(ratomFull)
					allRA.add(rkeys)
					la2raTbl[latomFull].add(rkeys)
			# second pass: freeze ratomSet; count all latomTypes associated with it
			for latomFull in la2raTbl.keys():
				rkeys = la2raTbl[latomFull]
				frozenRKeys = frozenset(rkeys)
									
				if frozenRKeys not in vdwpTbl:
					vdwpTbl[frozenRKeys] = {}
					
				# NB: selection of RASets based on particular latomFull, but frequencies aggregated based on latomTYPE
				latype = getAtomType(latomFull)
				if latype not in vdwpTbl[frozenRKeys]:
					vdwpTbl[frozenRKeys][latype] = 1
				else:
					vdwpTbl[frozenRKeys][latype] += 1
	
		print 'analRLIF: %s NRatoms=%d NVDWRASets=%d' % (exptName, len(allRA),len(vdwpTbl))

		## Accumulate frequencies of reference to vdwP across ligand atoms
		totCnt = {}
		allRASets = vdwpTbl.keys()
		# vdwPatchTbl: raSet -> latype -> freq
		for raset in allRASets:
			tot = sum([vdwpTbl[raset][latom] for latom in vdwpTbl[raset]])
			totCnt[raset] = tot
	
		# 2do: EXPT: test various clustering thresh for vdwPatch
		# simThresh=simThresh,
		# distThresh=distThresh
		
		# 2do ASAP: WEIGHT individual vdwP by their FREQUENCIES prior to clustering!
		
		vdwpSimTbl =  defaultdict(dict) # vdwp1 -> vdwp2 -> tsim
		nlofreq = 0
		for vdwp1 in allRASets:
			if totCnt[vdwp1] < config.minVDWPRef: 
				nlofreq += 1
				continue
			
			for vdwp2 in allRASets:
				if vdwp2 <= vdwp1:
					continue
	
				if totCnt[vdwp2] < config.minVDWPRef: 
					nlofreq += 1
					continue
			
				vdwpsim = jaccardSim(vdwp1,vdwp2)
				vdwpSimTbl[vdwp1][vdwp2] = vdwpsim
				vdwpSimTbl[vdwp2][vdwp1] = vdwpsim
				
		print 'analRLIF: %s NLoFreqVDWP=%d NFreqVDWP=%d' % (exptName,nlofreq,len(vdwpSimTbl))
	
		# vdwpClustTbl: cliqueIdx -> (ctrFrag, [cliqueFrags] )
		vdwpClustTbl = bldFragClusters(vdwpSimTbl)
	
		print 'analRLIF: %s NRASets=%d NClust=%d' % (exptName,len(allRASets),len(vdwpClustTbl))
		
		outf = VDWPDir + '%s_vdwpCtr.csv' % (exptName)
		outs = open(outf,'w')
		outs.write('Idx,CVDWPP,ResList,Tot,CtrVDWP,VDWPP,Freq,VDWP\n')
		for cidx in vdwpClustTbl:
			cvdwpl = list(vdwpClustTbl[cidx][0])
			cvdwpl.sort()
			cvdwps = str(cvdwpl)
			cvdwpp = vdwpPP(cvdwpl)
			resList = vdwpRes(cvdwpl)
			tot = 0
			for vdwp in vdwpClustTbl[cidx][1]:
				tot += totCnt[vdwp]
			for vdwp in vdwpClustTbl[cidx][1]:
				freq = totCnt[vdwp]
				vdwpl = list(vdwp)
				vdwpl.sort()
				vdwps = str(vdwpl)
				vdwpp = vdwpPP(vdwpl)
				outs.write('%s,"%s","%s",%d,"%s","%s",%d,"%s"\n' % (cidx,cvdwpp,resList,tot,cvdwps,vdwpp,freq,vdwps))
		outs.close()
		
		## FINALLY associate any ligand's atoms which have vdw interaction with any ratom in ratomSet
		##		associate with ratom CLUSTER CENTROID
		
		# build dict for all vdwRASets -> centroid RASet
		vdraCtrTbl = {}  # vdwRASet -> centroid RASet
		# vdwpClustTbl: cliqueIdx -> (ctrFrag, [cliqueFrags] )
		for cidx in vdwpClustTbl.keys():
			for vdwRASet in vdwpClustTbl[cidx][1]:
				vdraCtrTbl[vdwRASet] = vdwpClustTbl[cidx][0]
  
		## First pass thru interTbl (via la2vdwTbl): add vdwPatch RLIF
				
		# la2rTbl: ligIdx -> latomFull -> [(rchain,raa,ratom)]
		nmissvra = 0
		for il,ligIdx in enumerate(la2vdwTbl.keys()):
			for latomFull in la2vdwTbl[ligIdx]:
				rfeatSet = la2vdwTbl[ligIdx][latomFull]
				rfeatList = ['_'.join(rfeat) for rfeat in rfeatSet]
				vdwRASet = frozenset(rfeatList)
				if vdwRASet not in vdraCtrTbl:
					# print 'analRLIF: vdwRASet not found?!',ligIdx,latomFull,vdwRASet
					nmissvra += 1
					continue
				CtrvdwRASet = vdraCtrTbl[vdwRASet] 
	
				allRes = vdwpPP(CtrvdwRASet) # "chain+residue:[resAtoms];chain+residue:[resAtoms];..."
				# latype = getAtomType(latomFull)
				# 151127: drop ligand details from vdw-RLIF; seems to carry little additional info
				rlif = bldFeature3(allRes,'vdw','vdw','vdw','vdw')
								
				# 151127: include latomFull in rlif index
				rlifTbl[rlif].append( (ligIdx,latomFull) )									  
		
		print 'analRLIF: %s NvdwPatch = %d NMissvdwRAtomSet=%d' % (exptName,len(rlifTbl),nmissvra)
	
	## Second pass thru interTbl: add non-vdw RLIF
	for rfeat,ligTbl in interTbl.items():			
		# (rchain,raa,ratom) -> {itype - > {latom -> [ ligID ] } }
		(rchain,raa,ratom) = rfeat
		for itype,ligFeatTbl in ligTbl.items():
			
			# ASSUME vdw handled above
			if itype == 'vdw' and config.vdwLigandAtom:
				continue

			for latomFull,ligIDList in ligFeatTbl.items():					
				f = bldFeature3(rchain,raa,ratom,itype,latomFull)
				for ligIdx in ligIDList:
					# 151127: include latomFull in rlif index
					rlifTbl[f].append( (ligIdx, latomFull) )
	
	print 'analRLIF: %s NNonZinc=%d NRLIF=%d' % (exptName,config.NNonZincLig,len(rlifTbl))	   
	
	allRLIF = rlifTbl.keys()
	
	# Save RLIF stats
	allRLIF.sort(key= lambda f: len(rlifTbl[f]),reverse=True)	   
	itypeFreqTbl = defaultdict(int)
	itypeCummTbl = defaultdict(int)
	
	outf = RLIFDir + ('%s_rlifFreq.csv' % exptName)
	outs = open(outf,'w')
	outs.write('RLIF,Freq\n')
	for f in allRLIF:
		outs.write('"%s",%d\n' % (f,len(rlifTbl[f])))			
		fbits = feature2bits(f) # [chain,RAA,Ratom,IType,Latom]
		itype = fbits[3]
		itypeFreqTbl[itype] += 1
		itypeCummTbl[itype] += len(rlifTbl[f])		 
	outs.close()
	
	cntStr = ''
	for itype in InterTypes:
		if itype in itypeFreqTbl:
			cntStr += ' %s=%d/%d' % (itype,itypeFreqTbl[itype],itypeCummTbl[itype])
		else:
			cntStr += ' %s=0/0' % (itype)
	print 'analRLIF: IType Distrib %s: %s' % (exptName,cntStr) 

	# Save RLIF2Lig
	allRLIF.sort()
	
	outs = open(r2lf,'w')
	outs.write('RLIF,Ligand,LAtom\n')
	for f in allRLIF:
		for ligIdx,latomFull in rlifTbl[f]:
			zincid = ligIdx2zinc(ligIdx)
			outs.write('"%s",%s,%s\n' % (f,zincid,latomFull))

def loadVDWPCtrTbl(vdwpcf):
	''' vdwPatchSet -> vdwpCtr
	'''
	reader = csv.DictReader(open(vdwpcf))
	vdwPCtrTbl = defaultdict(list)
	for i,entry in enumerate(reader):
		# Idx,CVDWPP,ResList,Tot,CtrVDWP,VDWPP,Freq,VDWP
		vdwp = eval(entry['VDWP'])
		vdwpSet = frozenset(vdwp)
		vdwpCtr = eval(entry['CtrVDWP'])
		vdwPCtrTbl[vdwpSet] = vdwpCtr
		
	print 'loadVDWPCtrTbl: NVDWP=%d' % (len(vdwPCtrTbl))
	return vdwPCtrTbl
	
def loadRLIF2LigTbl(r2lf):
	'''return rlif2LigTbl: rlif -> [(ligIdx,latom)] from r2l file created by analRLIF()
	'''	  

	reader = csv.DictReader(open(r2lf))
	r2lTbl = defaultdict(list)
	for i,entry in enumerate(reader):
		# RLIF,Ligand,LAtom
		zincid = entry['Ligand']
		ligIdx = normLigand(zincid)
		rlif = entry['RLIF']
		latom = entry['LAtom']
		r2lTbl[rlif].append( (ligIdx,latom) )
		
	print 'loadRLIF2LigTbl: NRLIF=%d' % (len(r2lTbl))
	return r2lTbl
		
def bldLig2RLIFTbl(r2lf):
	'''return inverted lig2rlifTbl: ligIdx -> rlif -> latom from r2l file created by analRLIF()
	'''	  
	 
	reader = csv.DictReader(open(r2lf))
	l2rTbl = defaultdict(dict)
	for i,entry in enumerate(reader):
		# RLIF,Ligand,LAtom
		zincid = entry['Ligand']
		ligIdx = normLigand(zincid)
		rlif = entry['RLIF']
		latom = entry['LAtom']
		l2rTbl[ligIdx][rlif] = latom
		
	print 'bldlig2RLIFTbl: NLig=%d' % (len(l2rTbl))
	return l2rTbl

def bldLigAtom2vdw(interTbl):
	
	## Build inverted dict ligIdx -> latomFull -> [(rchain,raa,ratom)]
	la2rTbl = defaultdict(lambda: defaultdict(list))
	for rfeat,ligTbl in interTbl.items():			
		# (rchain,raa,ratom) -> {itype - > {latom -> [ ligID ] } }
		(rchain,raa,ratom) = rfeat
		for itype,ligFeatTbl in ligTbl.items():
			
			if itype != 'vdw':
				continue
			for latomFull,ligIDList in ligFeatTbl.items():
				for ligIdx in ligIDList:
					la2rTbl[ligIdx][latomFull].append(rfeat)
	return la2rTbl

def loadNonZinc(exptName):
	config.NonZincLigTbl = {} # ligand name -> nonzincID 
	config.NZIdx2LigTbl = {} # nonzincID -> ligand name
	config.NNonZincLig = 0
	nonzf = LowEDir + 'nonZincLig.csv'

	if os.path.exists(nonzf):
		# cf. FAAHA_expt()
		reader = csv.DictReader(open(nonzf))
		for i,entry in enumerate(reader):
			# Ligand,NZIdx
			config.NonZincLigTbl[ entry['Ligand'] ] = int(entry['NZIdx'])
			config.NZIdx2LigTbl [ int(entry['NZIdx']) ] =  entry['Ligand']
		config.NNonZincLig = len(config.NonZincLigTbl)
		print 'loadNonZinc: %d NonZincLig loaded from %s' % (config.NNonZincLig,nonzf)
	else:
		print 'loadNonZinc: no NonZincLig',exptName
	
def loadBestLig(exptName,lowef,maxLig=-1,normalizeLigand=True):
	'''read energy, batch from lowef
	return ligTbl: ligIdx -> [e,batch] and batch2ligTbl: batch -> [ligs in batch]
	if maxLig != -1, only maxLig lowest energy ligands included
	150629: normLigand() applied
	151028: preload NonZincLigTbl, ligTbl assumes it
	'''

	if normalizeLigand:
		loadNonZinc(exptName)
		
	ligTbl = {} # lig -> [e,batch]																							  

	reader = csv.DictReader(open(lowef))
	for i,entry in enumerate(reader):
		# cf FAAHA																									   
		#  fldNameLine = 'Batch,Ligand,E,Eff,Nvdw,Ninter\n'  
																	   
		# HACK: let batch tag along in ligTbl in first pass
		ligand = entry['Ligand']
		batchNo = int(entry['Batch'])
		
		if normalizeLigand:
			ligIdx = normLigand(ligand)													 
			ligTbl[ ligIdx ] = [ float( entry['E'] ), batchNo ]
		else:
			ligTbl[ ligand ] = [ float( entry['E'] ), batchNo ]

	batch2ligTbl = defaultdict(list) # batch -> [lig]																	   
	if maxLig==-1:
		for lig in ligTbl:
			batch2ligTbl[ ligTbl[lig][1] ].append(lig)
		return ligTbl, batch2ligTbl

	lowestLig = ligTbl.keys()
	# ASSUME: natural ordering correct for these negative energies														  
	lowestLig.sort(key=lambda k: ligTbl[k][0])

	# remove larger energy ligands
	for bigLig in lowestLig[maxLig:]:
		del ligTbl[bigLig]

	for lig in ligTbl:
		batch2ligTbl[ ligTbl[lig][1] ].append(lig)
			
	print 'loadBestLig: %s %d ligands loaded, %d batches' % (exptName,len(ligTbl),len(batch2ligTbl))
	return ligTbl, batch2ligTbl

def selectEliteDist(faahDir,exptKey,batchList,thresh,rabbleFrac,rqual,runType='ADV',dcrit='energy'):
	''' designed to be plug-compatible with loadBestLig()
	Create sample of ELITE (e<thresh) + Rabble=random non-elite of rabbleRatio size wrt/ ELITE
	rabbleFrac=thresh means ~ same number of each; rabbleFrac=1.0 means ALL non-elite included
	NB: requires exhaustive review of all ligands across crawls summ files
	NB: only ligands in l2rTbl considered (eg, if no vdwFeatures, wont be included; considered skipped)
	requires getting processed files out of tarball
	write eliteDist
	ligIdx -> [e,batch] and batch2ligTbl: batch -> [ligs in batch]
	
	'''

	exptNo,prot,recept,site,lib = exptKey															  
	exptName = bldExptStr(exptKey)
	
	loadNonZinc(exptName)

	allE = []
	noddE = 0
	nmissf = 0
	nrcd=0
	ndup=0
	ligTblE = {} # ligidx -> [batchNo,dval]
	
	for isf,batchNo in enumerate(batchList):
		summPath = faahDir+('summ/%s_summ_%07d.csv' % (runType,batchNo))
			
		try:
			inStr = open(summPath)
		except Exception,e:
			print "getThresh: can't open1 %s?!" % (summPath)
			nmissf += 1
			continue
		
		for il,line in enumerate(inStr.readlines()):
			if il == 0:
				# Expt,Recept,Ligand,E,Eff,Nvdw,Ninter
				# written by crawl_ADV.rptData_ADV()
				continue

			nrcd += 1
			try:
				flds =line[:-1].split(',')
				(expt,batch,ligRecept,ligand,e,eff,nvdv,ninter) = flds
					
				# 160109: this test is getting in the way of slight change in crawl output?
#				 if ligRecept != exptRecept:
#					 continue
				
				## 141219: Pass 1: selection criterion made a parameter
				
				if dcrit == 'energy':
					dval = float(e)
				elif dcrit == 'ligEff':
					dval = float(eff)
				elif dcrit == 'ligEnth':
					# Assume LigDOFTbl already loaded
					if ligand not in LigDOFTbl:
						print 'getThresh: lig missing from LigDOFTbl?!',ligand
						continue
					dof = LigDOFTbl[ligand]
					natom = float(e) / float(eff)
					natomi = int(natom)
					assert natomi-natom > 1e-2, 'huh?'
						
					dval = (float(e) + (dof * FE_coeff_tors)) / natom	
											
			except Exception, e:
				print "getThresh: bad line(1) %s %d %s?!" % (summPath,il,e)
				noddE += 1
				break # don't try to read other lines from inStr
			
			# 160122: Remove redundant ligands
			ligIdx = normLigand(ligand)
			if ligIdx in ligTblE:
				# NB: dups written by getThresh(initRun=True)
				ndup += 1
				continue

			ligTblE[ligIdx] = [batchNo,dval]

			# Pass1: need to make this thresh criterion-sensitive, too
			
			if (dcrit == 'energy' and abs(dval) > TooBigE) or \
				(dcrit == 'ligEff' and abs(dval) > TooBigLigEff) or \
				(dcrit == 'ligEnth' and abs(dval) > TooBigLigEnth):
				noddE += 1
				continue
								
			allE.append(dval)
			
		inStr.close() # eo-batch file

	nlig = len(ligTblE)
	print 'selectEliteDist: %s NLig=%d NDup=%d NOddE=%d' % (exptName,nlig,ndup,noddE)
		
	# 160103: what's this for?!
	# lowestLig = ligTbl.keys()
	# ASSUME: natural ordering correct for these negative energies														  
	# lowestLig.sort(key=lambda k: ligTbl[k][0])

	batch2ligTbl = defaultdict(list) # batch -> [ligand]  

	outDir = EliteDir
	eligf = EliteDir + exptName + '_eliteLig.csv'
	outs = open(eligf,'w')
	outs.write('Ligand,Batch,E,Elite\n')
	nlowe = 0
	nrabble = 0
	nother = 0
	nmissLig = 0
	RabbleMoatEnergyGap = 0.2 # 160523: leave gap between elite and rabble energies
	for il,ligIdx in enumerate(ligTblE):
		zincid = ligIdx2zinc(ligIdx)
		# 160714
		# Bad zincid, perhaps because batch is missing from process2?
		if zincid == 'Lig??':
			nmissLig += 1
			continue
			
		bnos = ligTblE[ligIdx][0]
		e = ligTblE[ligIdx][1]
		if e < thresh:	
			batch2ligTbl[ bnos ].append(zincid)
			nlowe += 1
			outs.write('%s,%s,%f,%d\n' %(zincid,bnos,e,1))
			
		# NB leave gap between elite and rabble energies
		# NB: rabbleFrac number of RABBLE
		elif e > (thresh + RabbleMoatEnergyGap) and random.random() < rabbleFrac:
			batch2ligTbl[ bnos ].append(zincid)  
			nrabble += 1
			outs.write('%s,%s,%f,%d\n' %(zincid,bnos,e,0))
		else:
			nother += 1
	outs.close()

	ligTblE = {} # can't del ligTbl?
	
	print 'selectEliteDist: %s Thresh=%f RabbleFrac=%f NLowE=%d NRabble=%d NOther=%d NMissLig=%d TotDistLig=%d' % \
		(exptName,thresh,rabbleFrac,nlowe,nrabble,nother,nmissLig,nlowe+nrabble)
	
	allBatch = batch2ligTbl.keys()
	allBatch.sort()

	outDir = EliteDir + exptName + '/'
	nligOut = 0
	for batch in allBatch:
		bnos = '%07d' % int(batch)
		ligList = batch2ligTbl[batch]
			
		# .../processed/Exp120/Results_x1ZTZ_prASw0c0/FAHV_x1ZTZ_prASw0c0_0550316_processed.tgz
		# tgzPath = ProcDir+'Exp%s/Results_%s/FAHV_%s_%s_processed.tgz' % (exptNo,rqual,rqual,bnos)
		
		# 160111: modified name of tgz, 
		#		 first check to see if uncompressed version exists
		# batchDir = ProcDir+'Exp%s/Results_%s/batch_%s/' % (exptNo,rqual,bnos)
		
		batchDir = ProcDir + '%s/%s/batch_%s/' % (exptNo,recept,bnos)
		
		if os.path.isdir(batchDir):
			getPDBQT.mgl_getADV_PDBQT(exptName,batchDir,bnos,ligList,outDir)
		else:
			tgzPath = ProcDir+'Exp%s/Results_%s/FAHV_%s_%s_processed.tgz' % (exptNo,rqual,rqual,bnos)
			getPDBQT.mgl_getADV_PDBQT_tgz(exptName,tgzPath,bnos,ligList,outDir)
			
		nligOut += len(ligList)

	print 'selectEliteDist: %s NEliteDistOut=%d' % (exptName,nligOut)

def provideEliteDist(exptName):
	'''return ligTbl, batch2ligTbl with contents of EliteDir
	NB: 
	'''

	inf = EliteDir + '%s_eliteLig.csv' % (exptName)
	reader = csv.DictReader(open(inf))
	ligTbl = {} # ligand -> e
	batch2ligTbl = defaultdict(list) # batch -> [lig]
	nelite = 0
	for i,entry in enumerate(reader):
		# Ligand,Batch,E,Elite
		e = float(entry['E'])
		batchNo = int(entry['Batch'])
		zincid = entry['Ligand']
		ligIdx = normLigand(zincid)
		elite = int(entry['Elite'])
		if elite:
			nelite += 1
		ligTbl[ligIdx] = (e,batchNo,elite)
		batch2ligTbl[batchNo].append(ligIdx)

	print 'provideEliteDist: %s NElite=%d NLig=%d' % (exptName,nelite,len(ligTbl))
	return ligTbl, batch2ligTbl
 
def provideEliteTrain(exptName):
	'''return set of elite ligands from EliteDir
	'''

	inf = EliteDir + '%s_eliteLig.csv' % (exptName)
	reader = csv.DictReader(open(inf))
	eliteSet =  set()
	nelite = 0
	nlig = 0
	for i,entry in enumerate(reader):
		# Ligand,Batch,E,Elite
		ligIdx = normLigand(entry['Ligand'])
		elite = int(entry['Elite'])
		nlig += 1
		if elite:
			eliteSet.add(ligIdx)
			nelite += 1

	print 'provideEliteActiveTrain: %s NElite=%d NLig=%d' % (exptName,nelite,nlig)
	return eliteSet
   
AARE = '([A-Z]+)([0-9]+)'
AAREPat = re.compile(AARE)

def splitRAA(raa):
	'split X999 -> pos, X'
	m = AAREPat.match(raa)
	try:
		(aaname,aapos) = m.groups()
	except:
		aapos = '999'
		aaname = raa
	iapos = int(aapos)
	return iapos, aaname

def bldRAA4(raa):
	'''999A: nicely formatted POSAA
	where AA is single character key in AADict
	'''

	iapos, aaname = splitRAA(raa)
	if len(aaname)==1:
		aa1 = aaname
	elif len(aaname)==3:
		if aaname in AADict:
			aa1 = AADict[aaname]
		else:
			aa1 = 'X'
	else:
		aa1 = 'X'

	# assert iapos < 1000, "iapos > 1000; '%03d' won't work"
	if iapos >= 1000:
		# print 'huh',raa, iapos
		# hack 141001!
		prefix = '%s_%s%s' % (chain,'999',aaname)
		return prefix
	
	aapos2 = '%03d' % iapos
	raa4 = '%s%s' % (aapos2,aa1)
	
	return raa4
		
def bldFeaturePrefix(chain,raa):
	'''C_999A: just chain + nicely formatted POSAA
	where AA is single character key in AADict
	'''

	iapos, aaname = splitRAA(raa)
	if len(aaname)==1:
		aa1 = aaname
	elif len(aaname)==3:
		if aaname in AADict:
			aa1 = AADict[aaname]
		else:
			aa1 = 'X'
	else:
		aa1 = 'X'

	# assert iapos < 1000, "iapos > 1000; '%03d' won't work"
	if iapos >= 1000:
		# print 'huh',raa, iapos
		# hack 141001!
		prefix = '%s_%s%s' % (chain,'999',aaname)
		return prefix
	
	aapos2 = '%03d' % iapos
	prefix = '%s_%s%s' % (chain,aapos2,aa1)
	
	return prefix
	
def getAtomType(ligatom):
	'''drop ligand atom index'''
	
	if len(ligatom)==0:
		latype = ''
	else:
		# NB: some atoms (eg, Cl) are > once character!
		lastAlpha = -1
		for i in range(len(ligatom)-1,-1,-1):
			if not ligatom[i].isdigit():
				lastAlpha = i+1
				break
		latype = ligatom[:lastAlpha]
	return latype

def getAtomIndex(ligatom):
	'''get only ligand atom index'''
	
	if len(ligatom)==0:
		laIdx = -1
	else:
		# NB: some atoms (eg, Cl) are > once character!
		lastAlpha = 0
		for i in range(len(ligatom)-1,-1,-1):
			if not ligatom[i].isdigit():
				lastAlpha = i+1
				break
		try:
			laIdx = int(ligatom[lastAlpha:])
		except:
			return -1
	return laIdx

def bldFeature3(chain,raa,ratom,itype,ligatom):
	'''v3: simplified logic for ALL itypyes features
	incorporte vdwPatches
	NB: unicode (from interact JSON) converted to str()
	'''

	latype = getAtomType(ligatom)

	# ASSUME raa already normalized in bldFocalInterTbl()
	f = '_'.join([chain,raa,str(ratom),str(itype),latype]) 

	return f

def feature2bits(f):
	bits = f.split('_')
	return bits

def analTrueEnergy(exptName,ligTbl,activeIdxSet):	
	'''analyze low energies with benefit of knowledge of activeIdxSet
	NB: dropped edecimal rounding?'''
	
	posEVec = []
	negEVec = []
	posEDistTbl = defaultdict(int)
	negEDistTbl = defaultdict(int)
	
	for ligIdx in ligTbl.keys():
					   
			e = float(ligTbl[ligIdx][0])
			if abs(e) > TooBigE:
				# print 'FAAHA: odd energy?!',isf,summF,il,e
				continue

			rEstr = rndEstr(e,ndig=1) # default edecimal
		   
			if ligIdx in activeIdxSet:
				posEVec.append(e)
				posEDistTbl[rEstr] += 1
			else:
				negEVec.append(e)
				negEDistTbl[rEstr] += 1

	posAvg,posSD = basicStats(posEVec)
	negAvg,negSD = basicStats(negEVec)
	print "analTrueEnergy: %s PosN=%d PosAvg=%f PosSD=%f NegN=%d NegAvg=%f NegSD=%f" % \
			(exptName,len(posEVec),posAvg,posSD,len(negEVec),negAvg,negSD)
			
	cummTrueMaxE, fracNegDropped = plot.plot2ExptE(exptName,posEDistTbl,negEDistTbl)
	return cummTrueMaxE, fracNegDropped

def lig2SpArff(exptName,ligTbl,origFragClustTbl,r2fTbl,activeIdxSet,arrfFile,\
			ethresh=None,useFragQual=True,useVDWP=True,splitRAtomAttrib=False,wgtLowE=1.0):
	'''create ARFF encoding of expts r2fTbl wrt/ RLIF+CentroidFrag qualification
	ASSUME r2fTbl already built, activeIdxSet identified
	use fragClustf to collapse fragments to cluster center
	output ARFF includes ligand names; filtered by weka
	'''
	
	# Build fragClustTbl: rlif -> frag -> fragCtr
	fragClustTbl = {}
	nloInfo = 0
	for rlif in origFragClustTbl:

		fragClustTbl[rlif] = {}
		for cidx in origFragClustTbl[rlif]:
			ctr = origFragClustTbl[rlif][cidx][0]
			for frag in origFragClustTbl[rlif][cidx][1]:
				fragClustTbl[rlif][frag] = ctr
		
		
	# allLig = ligCoordTbl.keys()
	# r2fTbl: RLIF -> frag -> [ (ligIdx,fragIdx,latomFull,fragInfo) ]
	nmissRLIF = 0
	l2rfTbl = defaultdict(list) # ligIdx -> [ (rlif,frag,fragIdx,latomFull) ]
	for rlif in r2fTbl.keys():
		if rlif not in fragClustTbl:
			nmissRLIF += 1
		for frag in r2fTbl[rlif]:
			for ligInfo in r2fTbl[rlif][frag]:
				(ligIdx,fragIdx,latomFull,fragInfo) = ligInfo
				# NB: dropping fragInfo because it is so big!
				l2rfTbl[ligIdx].append( (rlif,frag,fragIdx,latomFull) )
				
	print 'lig2SpArff: NLig = %d NMissFragClustRLIF=%d NLoInfoRLIF=%d' % (len(l2rfTbl),nmissRLIF,nloInfo)
	if ethresh != None:
		print 'lig2SpArff: ethresh=%f used' % (ethresh)

	allLig = l2rfTbl.keys()
	allLig.sort()
	
	## identify RAtom, RFQ features
	
	allRAtomsSet = set()  # allows enumeration when coding data

	# both ratomLigTbl and rfqTbl are dict for speedy ligand lookups
	ratomLigTbl = defaultdict(dict) # ra -> {lig: True}
	rfqTbl = defaultdict(dict) # ratom_frag_itype_latom -> {lig: True}
	nmissFrag = 0
	nmissrlif = 0
	nfndFrag = 0
	nhiELig = 0
	nlig = 0
	nmissRLIF2 = 0
	nloInfo2 = 0
	missFragTbl = defaultdict(int)
	for il,ligIdx in enumerate(allLig):
		
		if ethresh != None:
			e = ligTbl[ligIdx][0]
			if e > ethresh:
				nhiELig += 1
				continue
		nlig += 1
		allMatch = l2rfTbl[ligIdx]
		allMatch.sort(key = lambda m: m[0] ) # sort on RLIF
		
		for mi,match in enumerate(allMatch):
			rlif = match[0]

			# add full but un-fragment-qualified RLIF
			allRAtomsSet.add(rlif)
			ratomLigTbl[rlif][ligIdx] = True

			# 151114: Add (non-vdwP, non-stack) RLIF without ligatom in all cases!
			# 160422: and not stack interactions

			if rlif.find('vdw') == -1 and \
				rlif.find('tpi') == -1 and rlif.find('ppi') == -1:
				bits = rlif.split('_')
				ra = '_'.join(bits[:4])			
				allRAtomsSet.add(ra)
				ratomLigTbl[ra][ligIdx] = True

				# if splitRAtomAttrib, create sub-RA = chain+res, chain+res+ratom
				if splitRAtomAttrib:
					cres = '_'.join(bits[:2])
					cresa = '_'.join(bits[:3])
					ratomLigTbl[cres][ligIdx] = True
					allRAtomsSet.add(cres)
					ratomLigTbl[cresa][ligIdx] = True
					allRAtomsSet.add(cresa)


#			 if rlif not in rlifInfoTbl:
#				 nmissRLIF2 += 1
#				 continue
#			 rlifInfo = rlifInfoTbl[rlif]
#			 itype = bits[3]
#			 if rlifInfo < config.RLIFInfoThreshTbl[itype]:
#				 nloInfo2 += 1
#				 continue
					   
			# 151130: experiments suggest using only UNQUALIFIED vdwPatches is best balance
			if not useFragQual or rlif.find('vdw') != -1:
				continue  

			frag = match[1]
			
			if rlif not in fragClustTbl or frag not in fragClustTbl[rlif]:
				# print 'ligRFC2arff: missing from fragClustTbl?!',il,lig,mi,rlif,frag
				if rlif not in fragClustTbl:
					nmissrlif += 1
				nmissFrag += 1
				missFragTbl[rlif + config.RFQJoinChar + frag] += 1
				continue
			
			nfndFrag += 1
			cfrag = fragClustTbl[rlif][frag]				
			
			rfq = rlif + config.RFQJoinChar + cfrag
			rfqTbl[rfq][ligIdx] = True
	
	if ethresh != None:
		print 'lig2SpArff: ethresh=%f ==> nhiELig=%d NLig=%d' % (ethresh,nhiELig,nlig)
	print 'lig2SpArff: %s NRA=%d NRFQ = %d NFndFrag=%d NMissFrag=%d NMissRLIF=%d NLoInfo2=%d NUniqMiss=%d' % \
		(exptName,len(ratomLigTbl),len(rfqTbl),nfndFrag,nmissFrag,nmissRLIF2,nloInfo2,len(missFragTbl))
	
	missFragFile = ArffDir+('%s_missFrag.csv' % (exptName))
	outs = open(missFragFile,'w')
	allMiss = missFragTbl.keys()
	allMiss.sort(key=lambda k: missFragTbl[k],reverse=True)
	outs.write('Miss,NMiss\n')
	for miss in allMiss:
		outs.write('%s,%d\n' % (miss,missFragTbl[miss]))
	outs.close()

	atIdx = {}

	# 160227:  HACK to exclude VDWP
	# arrfFile = ArffDir+('%s_noVDWP_sp.arff' % (exptName))
	# arrfFile = ArffDir+('%s_sp.arff' % (exptName))
	
	outs = open(arrfFile,'w')
		
	outs.write('@relation %s\n' % (exptName))
	
	# http://weka.wikispaces.com/ARFF+(stable+version) 
	# 
	# the first string value is assigned index 0: this means that,
	# internally, this value is stored as a 0. When a SparseInstance is
	# written, string instances with internal value 0 are not output, so
	# their string value is lost (and when the arff file is read again, the
	# default value 0 is the index of a different string value, so the
	# attribute value appears to change). To get around this problem, add a
	# dummy string value at index 0 that is never used whenever you declare
	# string attributes that are likely to be used in SparseInstance objects
	# and saved as Sparse ARFF files.

	outs.write('@attribute "0000-dummyVar" string\n')
	atIdx['dummyVar'] = 0

	outs.write('@attribute "0001-ligand" string\n')
	atIdx['ligand'] = 1
	
	outs.write('@attribute "0002-e" numeric\n')
	atIdx['e'] = 2

	## RAtom features
	
	prevRAtom = None
	
	allRAtomsList = list(allRAtomsSet)
	allRAtomsList.sort()
	
	nlowFreqRLIF = 0
	nvdwp = 0
	begIdx = len(atIdx)
	allRAtomUsedList = []  # separate list, because some of rlif in allRAtomsList will be disqualified

	for ir,ra in enumerate(allRAtomsList):
		# possibily filter vdwPatches
		if not useVDWP and ra.find('_vdw') != -1:
			continue
			
		if len(ratomLigTbl[ra]) < config.MinLigPerAttrib:
			nlowFreqRLIF += 1
			# print 'dropRLIF',len(ratomLigTbl[ra]),ra
			continue
		
		allRAtomUsedList.append(ra)
		atIdx[ra] = begIdx
	
		# truncate vdw suffix from ARFF attribute names
		if ra.find('_vdw') != -1:
			ra = ra.replace('_vdw_vdw_vdw_vdw', '')
			nvdwp += 1
				
		# NB: quotes required on attribute names, since these can contain commas!
		outs.write('@attribute "%04d-%s" {0,1}\n' % (begIdx,ra) )

		begIdx += 1
		
	## RFQ features	  
	allRFQ = rfqTbl.keys()
	allRFQ.sort()

	## Cull lowFreq RFQ, RFQ based on lowInfo RLIF
	nrfqFeature = 0
	nlofreqRFQ = 0
	nloInfoRLIF = 0
	allRFQUsedList = []  # separate list, because some of rfq in allRFQ will be disqualified
	begIdx = len(atIdx)
	for ir, rfq in enumerate(allRFQ):
		bits = rfq.split(config.RFQJoinChar)
		assert len(bits)==2, 'Odd RFQ %d %d %s %s' % (ir, len(bits),bits, rfq)
		[rlif,frag] = bits

		# drop low frequency fragments
		if len(rfqTbl[rfq]) < config.FragMinLigFreq:
			nlofreqRFQ += 1
			# print 'dropRFQ',len(rfqTbl[rfq]),rfq
			continue
		
#		 if rlif not in hiRLIFTbl:
#			 nloInfoRLIF += 1
#			 continue
			
		nrfqFeature += 1
		allRFQUsedList.append(rfq)   
		# NB: quotes required on attribute names, since these can contain commas!	
		outs.write('@attribute "%04d-%s" {0,1}\n' % (begIdx,rfq) )
		atIdx[rfq] = begIdx
		begIdx += 1

	print 'lig2SpArff: %s NRAtom=%d NRFQ=%d NLoFreqRLIF=%d NLoFreqRFQ=%d NLoInfoRLIF=%d NvdwP=%d' % \
		(exptName,len(allRAtomsSet),len(allRFQUsedList),nlowFreqRLIF,nlofreqRFQ,nloInfoRLIF,nvdwp)
	
	outs.write('@attribute class {0,1}\n')
	atIdx['class'] = len(atIdx)
	
	outs.write('@data\n')
	
	nout = 0
	nactive = 0
	nhighE = 0
	nNoFeature = 0
	nbitsVec = []
	nhiELig = 0
	
	# lig2code will contain residual of ligands without any features
	lig2codeSet= set(ligTbl.keys())

	for il,ligIdx in enumerate(ligTbl):
#		 if il % 1000 == 0:
#			 print 'ligRFC2arff: writing data...',il
		zincid = ligIdx2zinc(ligIdx)	
		e = ligTbl[ligIdx][0]
		if ethresh != None:
			if e > ethresh:
				nhiELig += 1
				continue

		ligStr = '{'

		ligStr += '%d "%s",' % (atIdx['ligand'],ligIdx)
			
		ligStr += '%d %f,' % (atIdx['e'],e)
					
		if ligIdx not in l2rfTbl:
			nNoFeature += 1
			
			nbitsVec.append(0)
			
		lig2codeSet.discard(ligIdx)	  
		
		nout += 1

		nbitsSet = 0
		for ra in allRAtomUsedList:
			
			if ligIdx in ratomLigTbl[ra]:
				try:
					ligStr += '%d 1,' % (atIdx[ra])
				except:
					print 'huh'			
				nbitsSet += 1
		 
		for rfq in allRFQUsedList:
			if ligIdx in rfqTbl[rfq]:
				ligStr += '%d 1,' % (atIdx[rfq])
				nbitsSet += 1
				
		nbitsVec.append(nbitsSet)
		
		# ligIdx = normLigand(ligand)
		# active = ligIdx in activeIdxSet
		active = ligIdx in activeIdxSet
		if active:
			nactive += 1
			# NB: never weighting of actives; they're rare!
			ligStr += '%d %d}\n' % (atIdx['class'],active)
		else:
			if wgtLowE==1.0:
				ligStr += '%d %d}\n' % (atIdx['class'],active)
			else:
				ligStr += '%d %d},{%f}\n' % (atIdx['class'],active,wgtLowE)
				
		outs.write(ligStr)

	ndrop1 = len(lig2codeSet)
	
	outs.close()		

	print 'lig2SpArff: %s NHiELig=%d/%d NRFQ=%d NRAtom=%d NLow=%d NNoFeature=%d NDrop1=%d Nout=%d NActive=%d/%d' % \
		(exptName,nhiELig,len(ligTbl),len(rfqTbl),len(allRAtomsList),len(ligTbl),nNoFeature,ndrop1,nout,nactive,len(activeIdxSet))

	print 'lig2SpArff: %s NRFQFeature=%d' % (exptName,len(allRFQUsedList))
		
	avg,sd = basicStats(nbitsVec)
	print 'lig2SpArff: %s NBitsSet Avg=%f SD=%f' % (exptName,avg,sd)   

def exptName2bits(exptName):
	
	bits = exptName.split('_')
	if len(bits)==6:
		# eg, 119_PR_x3KFS_prASw0c0_AS_VM
		(exptNo,protein,recept,exptVar,site,lib) = bits
	elif len(bits)==5:
		# eg, 119_PR_x3I8WBprASw0c0_AS_VM
		(exptNo,protein,r_ev,site,lib) = bits
		recept=r_ev[:6]
		exptVar=r_ev[6:]
	elif len(bits)==8:
		# eg 48_IN_x3ZT1_A_IN_LEDGF_LEDGF_EN
		(exptNo,protein,recept,chain,prot2,site,site2,lib) = bits
		exptVar = '' # chain  # 2do: stick chain in exptVar?
	else:
		print 'exptName2bits: odd exptName?!',exptName
		return []
		
	return [exptNo,protein,site,recept,lib,exptVar]

def exptNameDiff(exptName1,exptName2):
	
	bits1 = exptName2bits(exptName1)
	bits2 = exptName2bits(exptName2)
			
	diffs = []
	for i,b in enumerate(bits1):
		if b != bits2[i]:
			diffs.append( (b, bits2[i]) )
		
	return diffs
														 
def getAtoms(mol):
	# http://openbabel.org/docs/dev/UseTheLibrary/PythonDoc.html#using-iterators
	
	# Note that OBMolTorsionIter returns atom IDs which are off by one. 
	# That is, you need to add one to each ID to get the correct ID. 
	
	# cf http://forums.openbabel.org/Read-pdb-files-in-C-td4657503.html
	# http://forums.openbabel.org/Unexpected-behavior-with-GetResidue-td4657736.html
	
	atomList = []

	tbl = ob.OBElementTable()

#	 atcntTbl = defaultdict(int)

	# print 'aidx,aindex,aid,atype,anum,aname1,aname2'

	for atom in ob.OBMolAtomIter(mol):
		aid = atom.GetId()
		atype = atom.GetType()
		r = atom.GetResidue()
		if r==None:
			aname = tbl.GetSymbol(atom.GetAtomicNum())
		else:
			aname = r.GetAtomID(atom).strip()

#		 aidx = atom.GetIdx()
#		 aindex = atom.GetIndex()
#		 anum = atom.GetAtomicNum()
#		 aname2 = tbl.GetSymbol(atom.GetAtomicNum())
			
#		 print aidx,aindex,aid,atype,anum,aname1,aname2
			
		atomList.append( (int(aid), aname, atype) )
	

	return atomList		

def mapLig2Features(ligIdx,fragments,ligMol,errs):

	obc = ob.OBConversion()
	obc.SetInFormat('smi')
				
	latoms = getAtoms(ligMol)
	
	# print 'LIGAND',zincid,latoms

	bindDict = defaultdict(dict) # fragIdx -> isoIdx -> [(a,b)]
	bindDict['latoms'] = latoms

	bindDict['nfrag'] = len(fragments)
	# print 'fragments',fragments
	
	for fi,frag in enumerate(fragments):
		fragMol = ob.OBMol()

		obc.ReadString(fragMol,frag)
		
		## TJO 151021
		fragMol.DeleteHydrogens()

		## TJO 151111
		# DeleteHydrogens() worked fine when using CompileMoleculeQuery
		# but using SMARTS matching, we need to string delete H atoms
		# Just fix up troublesome [nH]
		fragSmarts = frag.replace("[nH]", "n")
			
		fatoms = getAtoms(fragMol)
	
		# print 'FRAGMENT',fatoms
		bindDict[fi]['fragment'] = fragSmarts
		bindDict[fi]['fatoms'] = fatoms
			
		# # TJO 11/6/2015
		# new code block using smart pattern matching and mapping
		# (vs. CompileMoleculeQuery(fragMol), OBIsomorphismMapper)
		# intended to produce same data structures and mapping, but
		# with fewer (no!) missing mappings
		pat = ob.OBSmartsPattern()
		if  pat.Init(fragSmarts) and pat.Match(ligMol):
			if pat.NumMatches() == 0:
				zincid = ligIdx2zinc(ligIdx)
				errs.write("no isomorphs %s %d %s\n" % (zincid,fi,fragSmarts))
				bindDict[fi]['maps'] = None
			else:
				maps = list()
				for p in pat.GetUMapList():
					maps.append( [(a,b-1) for (a,b) in enumerate(p)] )
				bindDict[fi]['maps'] = maps
		else:
			zincid = ligIdx2zinc(ligIdx)
			# included in bldLig2frag error reporting
			# errs.write("error making mapper  %s %d %s\n" % (zincid,fi,frag))
			bindDict[fi]['maps'] = None
				
	return bindDict

def assignFrag(ligIdx,bindDict,errs,exptName):
	'''find consistent assignment of fragments to ligand
		across alternative isomorphic maps for each fragment
		~ branch&bound, beginning with least ambiguous fragments
		  (ie fewest maps), and minimizing the number of fragment
		  atoms not matching ligand atoms
		  
		NB: there may be unassigned fragments, despite their having
			come from the ligand:
			
		Subject:	 Re: recapping - PS
		Date:	 Wed, 25 Mar 2015 08:50:27 -0700
		From:	 TJ ODonnell <tjo@acm.org>		
		
		The match/map pays attention to stereochemistry.  The fragment stereo
		carbon [C@@H] does not match the (we know) corresponding non-stereo C
		in the input molecule.  Changing the fragment stereo atom does the
		trick.  Stereo smiles (or smarts) matching should always be done with
		care and only when you must have the exact stereo match.  This is one
		reason I remove stereochemistry in recap, since fragmentation can
		happen at stereocenters.
		'''

	frag2check = range(bindDict['nfrag'])
	# begin with LONGEST, least ambiguous fragments
	def mapLen(f):
		if bindDict[f]['maps'] == None:
			return 0
		else:
			return len(bindDict[f]['maps'])
			   
#	 frag2check.sort(key = (lambda f: len(bindDict[f]['maps'])))
	frag2check.sort(key = mapLen)

	availLAtoms = set([lidx for (lidx,lname,ltype) in bindDict['latoms'] ])
	fragBind = {}
	completeBind = True
	missedFrag = []
	for fi,f in enumerate(frag2check):
		if bindDict[f]['maps'] == None:
			continue
		
		bestMap = None
		maxFree = 0
		for im,amap in enumerate(bindDict[f]['maps']):
			ligBindings = [lidx for fidx,lidx in amap]
			
			nfree = sum( [1 for lidx in ligBindings if lidx in availLAtoms] )
			if nfree > maxFree:
				maxFree = nfree
				bestMap = amap

		if not bestMap:
#			 print '\n* incomplete assignment %s FragIdx=%d' % (zincid,f)
			completeBind = False
			missedFrag.append(f)
			nmiss = len(bindDict[f]['fatoms'])
			fragBind[f] = {'fragment': bindDict[f]['fragment'],
						   'mapIdx': None,
						   'nmiss': nmiss,
						   'useMap': None,
						   'flaNames': None,
						   'dropped': []}
			continue
		
		# capture parallel atomName labels associated with bestMap
		## ASSUME order maintained!
		
		## include ALL (including dropped) flaNames kept from bestMap
		flaNames = []
		for fidx,lidx in bestMap:
			(lidx,lname,ltype) = bindDict['latoms'][lidx]
			flaNames.append(lname)

		if maxFree == len(bestMap):
			ligBindings = [lidx for fidx,lidx in bestMap]
			dropped = []
			useMap = bestMap
		else:
			ligBindings = []
			dropped = []
			useMap = []
			for fidx,lidx in bestMap:
				if lidx in availLAtoms:
					ligBindings.append(lidx)
					useMap.append( (fidx,lidx) )
				else:
					# NB insert None into useMap when dropped
					useMap.append( (fidx,None) )
					dropped.append( (fidx,lidx) )
			
		availLAtoms = availLAtoms - set(ligBindings)
		mapFnd = True
		nmiss = len(bestMap)-maxFree
		
#		 laIdxList = [lidx for fidx,lidx in useMap]
#		 for (lidx,lname,ltype) in bindDict['latoms']:
#			 if lidx in laIdxList:
#				 flaNames.append(lname)
			
		fragBind[f] = {'fragment': bindDict[f]['fragment'],
					   'mapIdx': im,
					   'nmiss': nmiss,
					   'useMap': useMap,
					   'flaNames': flaNames,
					   'dropped': dropped}
		
	if not completeBind:
		errMsg = '"%s" "[ ' % (bindDict['canon'])
		for f in missedFrag:
			errMsg += '(%d, ""%s"")' % (f,bindDict[f]['fragment'])
		errMsg += ' ]"'
		zincid = ligIdx2zinc(ligIdx)
		errs.write('incomplete bind %s %s %s\n' % (exptName,zincid,errMsg))
	
	return fragBind

def breakMacroCycles(obmol, minRingSize=7,verbose=False):
	'''break macrocycles
	 minRingSize=7 can catch cyclic ureas
	 S. Forli, 9 Mar 16
	'''

	rings = obmol.GetSSSR()  
	# blue-obelisk:findSmallestSetOfSmallestRings
	# which implementation?
	# [BGdV04a] Berger, F. and Gritzmann, P. and De Vries, S.. Minimum cycle bases for network graphs, Algorithmica. no. 1, . 2004, pp. 51-62.
	# [FIG96] Figueras, J.. Ring Perception Using Breadth-First Search, J. Chem. Inf. Comput Sci.. vol. 36. 1996, pp. 986-991. 

	if verbose: print "breakMacroCycles: total rings found", len(rings)
	for r in rings:
		try:
			size = int(r.Size())
		except Exception,e:
			print 'breakMacroCycles: cant get rings?!',e
			return None
			
		if not size >= minRingSize:
			if verbose: print "breakMacroCycles: ring skipped (no macrocycle: %d-membered)" % size
			continue
		if verbose: print "breakMacroCycles: MACROCYCLE %d-membered" % size
		atoms = list(r._path)
		atoms.append( atoms[0] )
		if verbose: print "breakMacroCycles: atoms in ring", atoms
		for i in range(size):
			bondIdx = [ atoms[0+i], atoms[1+i] ]
			if verbose: print "breakMacroCycles: processing bond index", bondIdx
			bond = obmol.GetBond(bondIdx[0], bondIdx[1])
			if bond.IsSingle():
				if verbose: print "breakMacroCycles: break bond", bondIdx				
				batom = bond.GetBeginAtom()
				eatom = bond.GetEndAtom()
				obmol.DeleteBond(bond)
				obmol.AddHydrogens(batom)
				obmol.AddHydrogens(eatom)
				break
			
	return obmol


def bldDUDEMol2Idx(protein):
	'''create index for Mol2 files: ligand -> mol2File
	uses shell `find` to be robust across subdirectory (eg, active/, decoy/ subdirectories
	'''

	topMol2Dir = config.Mol2Dir + protein
	cmdList = ['find', topMol2Dir, '-name', '*.mol2']
	outstr = subprocess.check_output(cmdList)
	fileList = outstr.split()
	ndup = 0
	
	mol2IdxTbl = {} # ligand -> mol2File
	for f in fileList:
		# /Data/sharedData/coevol-HIV/WCG/DUDE/ligands_mol2/HIVRT/mol2/active/CHEMBL100056.mol2
		spos = f.rfind('/')
		ppos = f.rfind('.')
		ligand = f[spos+1:ppos]
		if ligand in mol2IdxTbl:
			print 'bldDUDEMol2Idx: dupe ligand?!',ligand,mol2IdxTbl[ligand],f
			ndup += 1
		mol2IdxTbl[ligand] = f

	print 'bldDUDEMol2Idx: NLig=%d NDup=%d' % (len(mol2IdxTbl),ndup)
	return mol2IdxTbl			

def applyRecap(amol,obc,pat):

	currRecap = recap.Recap(amol,4)
	for si,bondName in enumerate(currRecap.bondNames):
		pat.Init(currRecap.smarts[bondName])
		currRecap.apply(pat, si)
	currRecap.decide_multiples()
	currRecap.split()
	
	ligRecap = obc.WriteString(amol,True)
	# this returns both the canonSmiles string, but also PDBQT file name?!
	recapStr = ligRecap.split()[0]

	fragments = recapStr.split('.')
	return fragments


def bldLig2frag(ligTbl,exptName,exptLigDir,verbose=False):
	'''return ligFragTbl: ligIdx -> [ {fragment,mapIdx,nmiss,useMap,flaNames,dropped} ]
	doing RECAP fragmentation directly from PDBQT 
	vs. separate canon,RECAP steps as in _v1
	'''

# 	global SmartsPat
# 	SmartsPat = ob.OBSmartsPattern()
	
	MaxAtomMissThresh = 1
	
	# NB: dont use defaultdict; needs to be serializable/pickled!
	l2fTbl = {} # ligIdx -> [ {fidx,fragment,map} ]
	
	l2fFile = L2FDir + exptName + '_lig2Frag.csv'

	outs = open(l2fFile,'w')
	outs.write('Ligand,FragIdx,Fragment,Map\n')

	errFile = L2FDir + exptName + '_l2f_err.txt'
	errs = open(errFile,'w')
		
	if verbose:
		# initialize BindPPFile
		vpps = open(config.BindPPFile, 'w')
		vpps.close()
		
	obc = ob.OBConversion()
	obc.SetInAndOutFormats('pdbqt','can')
	obc.SetOptions("-i", obc.OUTOPTIONS) # produce smiles without isomeric or stereo information, TJO, 9 Apr 15
	
	obcMol2 = ob.OBConversion()
	obcMol2.SetInAndOutFormats('mol2','can')
	obcMol2.SetOptions("-i", obcMol2.OUTOPTIONS) # produce smiles without isomeric or stereo information, TJO, 9 Apr 15

	ngood=0
	npoor=0
	nslide=0
	nerr = 0
	neqMol2PDBQT = 0
	neqBreak = 0
	
	pat = ob.OBSmartsPattern()
	
	mcopener = macrocycleopener.MacrocycleOpener()
	
	for iz,ligIdx in enumerate(ligTbl.keys()):
			
		e,batch,elite = ligTbl[ligIdx]
			
		# zincid (vs ligIdx) used for error messages
		zincid = ligIdx2zinc(ligIdx)
				
		bno = int(batch)
		# NB: ProcDir fully qualified in dockFile()
		pdbqf = dockFile(exptName,bno,ligIdx,exptLigDir)

		if pdbqf == None:				
			errs.write('missing pdbqt file: %d,%s\n' % (bno,ligIdx))
			dls = open(config.DropLigFile, 'a')
			dls.write('%s,"bldLig2Frag: missing pdbqt file: %d,%s"\n' % (zincid,bno,ligIdx))
			dls.close()
			nerr += 1
			import pdb; pdb.set_trace()
			continue
		
		ligMol = ob.OBMol()
		obc.ReadFile(ligMol,pdbqf)

		# NB RECAP works directly from ligmol	
		# don't need to build canonical string
		# except to include it in bindDict
		
		ligSmilesC = obc.WriteString(ligMol)
		# WriteString returns both the canonSmiles string, but also PDBQT file name?!
		canon = ligSmilesC.split()[0]

		# RFQJoinChar used to join fragment to RLIF; better not be in fragment's canonSMILES!
		assert config.RFQJoinChar not in canon, 'RFQJoinChar found in canonSMILES?!'
		
		# NB: need to keep track whether pdbqt or Mol2 format used to build ligMol
		obc2use = obc
		
		## use Mol2 version of ligand if it exists
		if config.CheckMol2Lig:
			ligMol2 = ob.OBMol()
			ligand = ligIdx2zinc(ligIdx)
			# ASSUME: config.mol2Idx built
			if ligand in config.mol2Idx:
				mol2f = config.mol2Idx[ligand]
				
				if os.path.exists(mol2f):
					obcMol2.ReadFile(ligMol2,mol2f)
					
					try:
						mol2ligSmilesC = obcMol2.WriteString(ligMol2)					
						mol2canon = mol2ligSmilesC.split()[0]
						if canon != mol2canon:
							# if verbose: print 'bldLig2Frag: pdbqt != mol2: %s,%d,"%s","%s"' % (ligand,neqMol2PDBQT,canon,mol2canon)
							# NB:  ligand from pdbqt file dropped in favor of mol2 version
							ligMol = ligMol2
							canon = mol2canon
							obc2use = obcMol2
						else:
							neqMol2PDBQT += 1
	
					except Exception, e:
						print 'bldLig2Frag: cant form canonSmiles for mol2 file?!',e,mol2f
				else:
					print 'bldLig2Frag: cant read mol2 file?!',mol2f
			else:
				print 'bldLig2Frag: ligand missing for ligidx?!',ligand

							
		# make copy in order to alter mol
		amol = ob.OBMol(ligMol)

		aligSmilesC = obc2use.WriteString(amol)
		acanon = aligSmilesC.split()[0]
		if acanon != canon:
			print 'bldLig2Frag: different canon',iz,zincid,canon,acanon
# 
# 		# HACK need pybel Mol (:
# 		pamol = pybel.Molecule(amol)
# 		pomol = pybel.Molecule(ligMol)
# 		afp = pamol.calcfp()
# 		ofp = pomol.calcfp()
# 		if afp != ofp:
# 			print 'different FP',zincid

		fragments = applyRecap(amol,obc2use,pat)		
					
		if config.BreakMacroCycles:
			mcopener.processMol(ligMol)
			
			bligSmilesC = obc2use.WriteString(ligMol)
			bcanon = bligSmilesC.split()[0]

# 			if verbose:
# 				if canon != bcanon:
# 					print 'bldLig2Frag: macroCycleBreak modified canon',zincid,canon,bcanon
			
			bfragments = applyRecap(ligMol,obc2use,pat)
				
# 			if verbose:
# 				fragSet = set(fragments)
# 				bfragSet = set(bfragments)
# 				if fragSet == bfragSet:
# 					neqBreak += 1
# 				else:
# 					symDiff = fragSet.symmetric_difference(bfragSet)
# 					print 'bldLig2Frag: bfrags differ: %s,%d,"%s","%s","%s"' % (zincid,neqBreak,list(symDiff),fragments,bfragments)
		
			fragments = bfragments
			
		###################################
		
		bindDict = mapLig2Features(ligIdx,fragments,ligMol,errs)

		###################################
		
		# NB: augment bindDict with canon here vs. in mapLig2Features, 
		# since obc created here?
		
		bindDict['canon'] = canon
		
		frag2check = range(bindDict['nfrag'])
		badFragMap = [fi for fi in frag2check if bindDict[fi]['maps'] == None]
		if len(badFragMap) >0:
			ligCanon = bindDict['canon']
			errMsg = '"%s" [ ' % bindDict['canon']
			for fi in badFragMap:
				errMsg += '(%d, "%s")' % (fi,bindDict[fi]['fragment'])
			errMsg += ' ]'
			errs.write('Missing maps: %s,%s\n' % (zincid,errMsg))
			dls = open(config.DropLigFile, 'a')
			dls.write('%s,bldLig2Frag: Missing maps: %s\n' % (zincid,errMsg))
			dls.close()
			npoor += 1
			continue
		
		###################################
		
		fragAssmt = assignFrag(ligIdx,bindDict,errs,exptName)

		###################################
		
		totMiss = 0
		missList = []
		for f in fragAssmt:
			totMiss += fragAssmt[f]['nmiss']
			missList.append(fragAssmt[f]['nmiss'])
			
		if verbose:
			vpps = open(config.BindPPFile, 'a')
			if totMiss > 0:
				vpps.write('* %s NMiss=%d %s\n' % (zincid,totMiss,missList))
			ppBindDict(zincid,bindDict,vpps)
			vpps.close()

		if totMiss > MaxAtomMissThresh:
			# print 'bldLig2Frag: poor ligand: iz=%d ligIdx=%d zincid=%s totMiss=%d' % (iz,ligIdx,zincid,totMiss)
			dls = open(config.DropLigFile, 'a')
			dls.write('%s,bldLig2Frag: poor ligand totMiss=%d\n' % (zincid,totMiss))
			dls.close()
			npoor += 1
			continue
		
		ngood += 1
				
		ligPost = 0
		for fi,f in enumerate(fragAssmt):
			outs.write('%s,%d,%s,"%s"\n' % (zincid,fi,fragAssmt[f]['fragment'],fragAssmt[f]['useMap']))

			# NB: allow a single missing atom to slide!?
			if fragAssmt[f]['nmiss'] <= MaxAtomMissThresh:
				if fragAssmt[f]['nmiss'] > 0:
					nslide += 1
				 
				ligPost += 1
				# NB: not using defaultdict, because l2fTbl needs to be serializable/pickled!
				if ligIdx in l2fTbl:
					l2fTbl[ligIdx].append( fragAssmt[f].copy() )
				else:
					l2fTbl[ligIdx] = [ fragAssmt[f].copy() ]
			
		if ligPost==0:
			print 'bldLig2Frag: unposted lig?: iz=%d ligIdx=%d zincid=%s totMiss=%d' % (iz,ligIdx,zincid,totMiss)

			dls = open(config.DropLigFile, 'a')
			dls.write('%s,bldLig2Frag: unposted lig totMiss=%d\n' % (zincid,totMiss))
			dls.close()
			
		if verbose and (iz % 1000 == 0):
			print 'bldLig2Frag: nlig=%d nerr=%d npoor=%d' % (iz,nerr,npoor)
					
			
							   
	print 'bldLig2frag: NLig=%d NMol2EqPDBQT=%d NEQBreak=%d NGood=%d NErr=%d NPoor=%d NSlide=%d' % \
		(len(ligTbl),neqMol2PDBQT,neqBreak,ngood,nerr,npoor,nslide)
		
	outs.close()
	errs.close()
	
	return l2fTbl

def fragDist(l2fTbl):
	'''returns frequency distribution of ligands across fragments in l2fTbl
	'''
	
	freqTbl = defaultdict(int)

	# ligFragTbl: ligIdx -> {fragment,mapIdx,nmiss,useMap,flaNames,dropped}
	for ligIdx in l2fTbl.keys():
		for fragInfo in  l2fTbl[ligIdx]:
			frag = fragInfo['fragment']
			freqTbl[frag] += 1
	
	return freqTbl

def rptFragDist(fragDistTbl,outf):
	freqItems = freqHist(fragDistTbl)
	
	outs = open(outf,'w')
	outs.write('Fragment,F\n')
	for frag,freq in freqItems:
		outs.write('%s,%d\n' % (frag,freq))
	outs.close()
   
def ppBindDict(zincid,bindDict,outs):
	outs.write('* %s NFrag=%d NLigAtom=%d\n' % (zincid,bindDict['nfrag'],len(bindDict['latoms'])))
	for f in range(bindDict['nfrag']):
		info = bindDict[f]
		if info['maps'] == None:
			mapLen = 0
		else:
			mapLen =  len(info['maps'])

		outs.write('Frag %d %s NFragAtom=%d NMaps=%d\n' % \
				   (f,info['fragment'],len(info['fatoms']), mapLen))

	outline = '   |   |'
	for latom in bindDict['latoms']:
		outline += '%03d|' % (int(latom[0]))
	outs.write(outline+'\n')
	outline = ' F | I |'
	for latom in bindDict['latoms']:
		outline += '%3s|' % (latom[1])
	outs.write(outline+'\n')

	lblChars = string.hexdigits

	for f in range(bindDict['nfrag']):
		info = bindDict[f]
		if info['maps'] == None:
			outs.write('%3d|  No map(:\n' % (f))
			continue
	
		for im,amap in enumerate(info['maps']):
			outArr = defaultdict()		 
			outStr = '%3d|%3d|' % (f,im)
			for fidx,lidx in amap:
				if fidx<len(lblChars):
					outArr[lidx] = ' %s ' % (lblChars[fidx])
				elif fidx<2*len(lblChars):
					outArr[lidx] = ' %s!' % (lblChars[fidx-len(lblChars)])
				else:
					outArr[lidx] = ' !!' 

			for li in range(len(bindDict['latoms'])):
				if li in outArr:
					outStr += '%s|' % (outArr[li])
				else:
					outStr += '   |'
			outs.write(outStr+'\n')
	  
def loadActives(actives_file):
	actives = []
	fs=open(actives_file,'r')
	for active in fs.readlines():
		actives.append(active.strip())
	fs.close()
	return actives

def dockFile(exptName,bno,ligIdx,exptLigDir):
	# .../Dock/1_RT_x2ZD1_RT_NNRTI_NNRTInADJ_NNRTI_DD/0000010_ZINC06556034_0_out_Vina_VS.pdbqt

	ligand = ligIdx2zinc(ligIdx)
				
	# 160523: update for modern FAAH
	# cf getPDBQT.mgl_getADV_PDBQT()
	# 	newfname = '%s_%s_Vina_VS.pdbqt' % (bnos,ligand)
	
	pdbqf =  exptLigDir + '%s/%07d_%s_Vina_VS.pdbqt' % (exptName,bno,ligand)
	
# 	else:
# 		fname = ('%07d_' % bno) + ligand + '_Vina_VS.pdbqt'	 
# 		pdbqf =  exptLigDir + fname
	
	if os.path.exists(pdbqf):
		return pdbqf
	else:
		print 'dockFile: pdbqf missing?!',pdbqf
		import pdb; pdb.set_trace()		
		return None

def rptR2FTbl(r2fTbl,outf):
	
	allRLIF = r2fTbl.keys()  # RLIF -> frag -> [ (zincid,fragIdx) ]
	allRLIF.sort()
	outs = open(outf,'w')
	outs.write('RLIF,Fragment,NLig,"[(ligIdx,FragIdx)]"\n')
	for rlif in allRLIF:
		for frag in r2fTbl[rlif]:
			outs.write('%s,%s,%d,"%s"\n' % (rlif,frag,len(r2fTbl[rlif][frag]),r2fTbl[rlif][frag]))
	outs.close()
	
def bldRLIF2frag(exptName,ligTbl,l2fTbl,l2rTbl,errf,dockFileTbl=None,verbose=False):
	'''Aggregate RECAP fragments associated with RLIF
	RLIF -> frag ->  [ (ligIdx,fragIdx,map,latomFull,fragInfo) ]
	
	dockFileTbl added to make combineLigLibTarget() work
	'''

	nmissFrag = 0
	nmissflig = 0
	nmissrlig = 0
	nmissinter = 0
	nflig = 0
	r2fTbl = {} # RLIF -> frag -> [ (ligIdx,fragIdx,latomFull) ] 

	ligList = ligTbl.keys()
	ligList.sort()
	
	nmissR2F = 0
	
	errs = open(errf,'w')
	nerr = 0
	nHCignore = 0
	nstack = 0
	
	for iz,ligIdx in enumerate(ligList):
		zincid = ligIdx2zinc(ligIdx)
		
		if ligIdx not in l2fTbl:
			# print 'bldRLIF2frag: ligIdx missing from fragments?!',zincid
			nmissflig += 1
			continue

		if ligIdx not in l2rTbl:
			# print 'bldRLIF2frag: ligIdx missing from r2f?!',zincid
			nmissrlig += 1
			continue
		
		nflig += 1
		
		ligFragOrig = l2fTbl[ligIdx]
		ligFrags = []
		# remove cruff from fragment list
		for fi,fragInfo in enumerate(ligFragOrig):
			if fragInfo == None:
				errs.write('no fragInfo %s %s\n' % (zincid,fi))
				nerr += 1
				continue
			if fragInfo['flaNames'] == None:
				errs.write('no flaNames %s %s %s\n' % (zincid,fi,fragInfo))
				nerr += 1
				continue
			ligFrags.append(fragInfo)

		e,batch,elite = ligTbl[ligIdx]

		anyLigMatch = False
		# only want to report 'latomFull not found' error once, not for every RLIF
		laFullPrevNotFound = False
		for rlif in l2rTbl[ligIdx]:
			
			[rchain,raa,ratom,itype,latype] = feature2bits(rlif)
									
			latomFull = l2rTbl[ligIdx][rlif]
			latomType = getAtomType(latomFull)
			
			if itype=='vdw' and latomType in config.VDWExcludedLigAtoms:
				nHCignore += 1
				continue
						
			# 160422: include stacking itypes
# 			if itype=='tpi' or itype=='ppi':
# 				nstack += 1
# 				continue

			if itype=='ppi' or itype=='tpi':
				# NB: no fragment info for stacking interactions
				if rlif not in r2fTbl:
					r2fTbl[rlif] = {}
				if 'stack' not in r2fTbl[rlif]:
					r2fTbl[rlif]['stack'] = []
					
				r2fTbl[rlif]['stack'].append( (ligIdx,-1,latomFull,{}) )
				fragFnd = True
				anyLigMatch = True
				nstack += 1
				continue
				
			fragFnd = False
			for fi,fragInfo in enumerate(ligFrags):
				
				if latomFull in fragInfo['flaNames']:
					## FOUND: (full ligand atom name from docking found in 
					##	fragInfo['flaNames'] associated with ligand's fragments by bldLig2Frag-> mapLig2Features -> getAtoms() ->
					##		OBRESIDUE.GetAtomID(atom).strip() OR OBElementTable.GetSymbol(atom.GetAtomicNum()) 
					## all fragInfo from l2fTbl carried forward
					frag = fragInfo['fragment']
					if rlif not in r2fTbl:
						r2fTbl[rlif] = {}
					if frag not in r2fTbl[rlif]:
						r2fTbl[rlif][frag] = []
					   
					r2fTbl[rlif][frag].append( (ligIdx,fi,latomFull,fragInfo) )
					fragFnd = True
					anyLigMatch = True
					break
				## 
				
			if not fragFnd:
				# only want to report 'latomFull not found' error once, not for every RLIF
				if laFullPrevNotFound:
					continue
				laFullPrevNotFound = True
				errs.write('latomFull not found %s %s %s\n' % (zincid,latomFull,itype))
				nerr += 1
				nmissFrag += 1
				continue
		
		if not anyLigMatch:
			errs.write('No RLIF-Frag interactions %s\n' % (zincid))
			nmissR2F += 1
			dls = open(config.DropLigFile, 'a')
			dls.write('%s,bldRLIF2frag: No RLIF-Frag interactions\n' % (zincid))
			dls.close()
			
		if verbose and iz % 1000 == 0:
			print 'bldRLIF2frag: iz=%d nunboundLigAtom=%d' % (iz,nmissFrag)
			
	## eo-ligList
	errs.close()
	
	print 'bldRLIF2frag: NLig=%d NMissFragLig=%d NMissRLIFLig=%d NMissFrag=%d NStack=%d NHCIgnore=%d NErr=%d NMissR2F=%d' % \
		(nflig,nmissflig,nmissrlig,nmissFrag,nstack,nHCignore,nerr,nmissR2F)
		
	return r2fTbl

def analTargetR2FC(exptName,r2fTbl,r2fcTbl,outf,verbose=True):
	'''creates RLIF2FragmentCenter maps for active targets
	with respect to existing config.UseR2FC = r2fcTbl fragment clusters 
	'''
	
	# r2fTbl: RLIF -> frag -> [ (ligIdx,fragIdx,latomFull,fragInfo) ]
	outs = open(outf,'w')
	outs.write('RLIF,FragIdx,Sim,Frag,CFrag\n')
	for rlif in r2fTbl.keys():
		if len(r2fTbl[rlif]) == 0:
			# NB: FragIdx<0 ==> rlif not in R2FC table, or without fragments
			outs.write('"%s",-1\n' % (rlif))

		for tfrag in r2fTbl[rlif].keys():
			# find closest fragCtr within r2fcTbl[rlif] clusters
			
			if rlif not in r2fcTbl or len(r2fcTbl[rlif])==0:
				print 'bldTargetR2FC: rlif not in R2FC or empty?!',rlif
				# NB: FragIdx<0 ==> rlif not in R2FC table, or without fragments
				outs.write('"%s",-2\n' % (rlif))
				continue

			closeFrag = findCloseFrag(tfrag,r2fcTbl[rlif],findSib=True)
			
			mol1 = pybel.readstring('can',tfrag)
			fp1 = mol1.calcfp()
			
			sim = 0. 
			closeIdx = None
			# r2fcTbl: rlif -> cliqueIdx -> (ctrFrag, [cliqueFrags] )
			
			for cidx in r2fcTbl[rlif]:
				
				cfrag = r2fcTbl[rlif][cidx][0]
				
				mol2 = pybel.readstring('can',cfrag)
				fp2 = mol2.calcfp()
				
				## uses OB's to compute Tanimoto similarity ('|') betweenn fp1 and fp2
				tsim = fp1 | fp2
				
				if tsim > sim:
					sim = tsim
					closeIdx = cidx
			
			# NB: FragIdx= 0 ==> R2F target - R2FC centroid info
			cfrag = closeFrag[0][0]
			sim = closeFrag[0][1]
			outs.write('"%s",0,%f,%s,%s\n' % (rlif,sim,tfrag,cfrag))
			
			for ic,sibfrag in enumerate(closeFrag[1:]):
				
				mol2 = pybel.readstring('can',sibfrag)
				fp2 = mol2.calcfp()
				
				## uses OB's to compute Tanimoto similarity ('|') betweenn fp1 and fp2
				tsim = fp1 | fp2

				# NB: FragIdx>0 ==> R2FC centroid - cluster members info
				outs.write('"%s",%d,%f,%s,%s\n' % (rlif,ic+1,tsim,tfrag,sibfrag))
				
	outs.close()

def findCloseFrag(findFrag,clustTbl,findSib=True):
	'''assuming clustTbl is r2fcTbl[rlif]
	returns [ (closeFrag, sim), ...] most similar to findFrag
	first element of list is centroid; if findSib, cluster siblings returned also
	'''
	mol1 = pybel.readstring('can',findFrag)
	fp1 = mol1.calcfp()
	
	sim = 0. 
	closeIdx = None
	
	allClust = clustTbl.keys()
	allClust.sort()
	for cidx in allClust:
		
		cfrag = clustTbl[cidx][0]
		
		mol2 = pybel.readstring('can',cfrag)
		fp2 = mol2.calcfp()
		
		## uses OB's to compute Tanimoto similarity ('|') betweenn fp1 and fp2
		tsim = fp1 | fp2
		
		if tsim > sim:
			sim = tsim
			closeIdx = cidx
	
	# NB: FragIdx= 0 ==> R2F target - R2FC centroid info
	if closeIdx == None:
		return []
	
	closeFrag = [ (clustTbl[closeIdx][0], sim) ]
	
	if findSib:
		ctr = clustTbl[closeIdx][0]
		for ic,cfrag in enumerate(clustTbl[closeIdx][1]):
			# NB: don't include centroid twice
			if cfrag == ctr:
				continue
			
			mol2 = pybel.readstring('can',cfrag)
			fp2 = mol2.calcfp()
			
			## uses OB's to compute Tanimoto similarity ('|') betweenn fp1 and fp2
			tsim = fp1 | fp2
	
			closeFrag.append( (cfrag,tsim) )
			
	return closeFrag

	
def bldR2FCtr(exptName,r2fTbl,rlifInfoTbl,outf,verbose=False):
	''' returns: r2fcTbl: rlif -> cliqueIdx -> (ctrFrag, [cliqueFrags] )
	with rlif-specific fragment clustering
	symmetric entries both kept for bldFragClusters()
	outputs fragClustf with all RLIFs clusters
	'''

	# cf. cluster tuning, 150714
	simThresh = 0.7
	distThresh=0.5
	
	def getSeqPos(rlif):
		# cf bldRLIF2frag()
		# rlif = '_'.join([chain,res,ratom,itype,latom])
		bits = rlif.split('_')
		res = bits[1]
		pos = int(res[:3])
		return pos 
	
	outs = open(outf,'w')
	outs.write('Expt,RLIF,NFrag,NLig,NSimFrag,NClust,NSec,Notes\n')
	outs.close()  # appended below
	
	allRLIF = r2fTbl.keys()
	allRLIF.sort()
	r2fcTbl = {}  # rlif -> cliqueIdx -> (ctrFrag, [cliqueFrags] )
	nemptySim=0
	
	nloInfo = 0  
	nmissRLIF = 0  
	# r2fTbl: RLIF -> frag -> [ (ligIdx,fragIdx,latomFull,fragInfo) ] 
	for ir,rlif in enumerate(allRLIF):
		bits = rlif.split('_')
		itype = bits[3]

		if rlifInfoTbl != None and rlif not in rlifInfoTbl:
			nmissRLIF += 1
			continue
			rlifInfo = rlifInfoTbl[rlif]

			if rlifInfo < config.RLIFInfoThreshTbl[itype]:
				nloInfo += 1
				continue
			
		allFrag = r2fTbl[rlif].keys()
		
		allFrag.sort()
		nlowFreq = 0
		ndropLig = 0

		begTime = datetime.datetime.now()
		
		fragSimTbl =  defaultdict(dict) # frag1 -> frag2 -> tsim
		totLig = 0
		for frag1 in allFrag:

			nfragLig = len(r2fTbl[rlif][frag1])
			totLig += nfragLig
			if nfragLig < config.FragMinLigFreq:
				nlowFreq += 1
				ndropLig += nfragLig
				continue
			
			if frag1 != 'stack':
				mol1 = pybel.readstring('can',frag1)
				fp1 = mol1.calcfp()
			
			for frag2 in allFrag:
				if frag2 <= frag1:
					continue
				
				nfragLig = len(r2fTbl[rlif][frag2])
				if nfragLig < config.FragMinLigFreq:
					continue

				if frag2 == 'stack':
					tsim = 1.0
				else:
					mol2 = pybel.readstring('can',frag2)
					fp2 = mol2.calcfp()
					
					## uses OB's to compute Tanimoto similarity ('|') betweenn fp1 and fp2
					tsim = fp1 | fp2
				
				fragSimTbl[frag1][frag2] = tsim
				# NB: symmetric entries both kept for bldFragClusters()
				fragSimTbl[frag2][frag1] = tsim
			 
		if len(fragSimTbl)==0:
			r2fcTbl[rlif] = {}
			outs = open(outf,'a')
			outs.write('%s,"%s",%d,%d,0,0,0,emptySim ndropLig=%d\n' % (exptName,rlif,totLig,len(allFrag),ndropLig))
			outs.close()
			nemptySim += 1
			continue 
			
		rlifClustTbl = bldFragClusters(fragSimTbl,distThresh=distThresh)
		# rlifClustTbl: cliqueIdx -> (ctrFrag, [cliqueFrags] )
		
		elapTime = datetime.datetime.now() - begTime
		
		outs = open(outf,'a')
		outs.write('%s,"%s",%d,%d,%d,%d,%s\n' % (exptName,rlif,totLig,len(allFrag),len(fragSimTbl),len(rlifClustTbl),elapTime.total_seconds()))
		outs.close()

		r2fcTbl[rlif] = rlifClustTbl
				
	# outs.close()
	
	print 'bldR2FCtr: %s NRLIF=%d NLoInfoRLIF=%d NEmptySim=%d' % \
		(exptName,len(r2fcTbl),nloInfo,nemptySim)
	
	return r2fcTbl

def bldFragClusters(fragSimTbl,distThresh=0.5):
	'''cluster fragments using fastcluster and DISTANCE flattening
	built from (RLIF-specific) fragSimTbl
	returns fragClustTbl: cliqueIdx -> (ctrFrag, [cliqueFrags] )
	NB: default values simThresh=0.7,distThresh=0.5; cf. score_bldClust_dist_150714.log
	'''

	def getCenter(fragList):
		'''identify ~centroid of fragment list
		calculate cummSim for each frag1 wrt/ all other frag2
		'''
		
		if len(fragList)==1:
			return fragList[0]
		
		fragList.sort()
		maxSim = 0
		bestFrag = None
		nsim = 0
		fragList.sort()
		tot2 = 0
		for frag1 in fragList:
			tot = 0
			for frag2 in fragList:
				if frag2 == frag1:
					continue
				nsim += 1
				sim = fragSimTbl[frag1][frag2]
				tot += sim
			if tot > maxSim:
				maxSim = tot
				bestFrag = frag1
			tot2 += tot
			
		if bestFrag==None:
			print 'getCenter: all sim==0.?!'
			# arbitrarily pick first frag
			bestFrag = fragList[0]
				
#		 print 'getCenter: NFrag=%d NSim=%d bestAvgSim=%f AvgSim=%f' % \
#			 (nfrag,nsim,maxSim/(nfrag-1),tot2/nsim)

		return bestFrag

	# http://danifold.net/fastcluster.html
	# The argument X is either a compressed distance matrix or a collection of 
	# N observation vectors in D dimensions as an (N x D) array. 
	# Apart from the argument preserve_input, the methods have the same input 
	# and output as the functions of the same name in the package scipy.cluster.hierarchy. 
	# Therefore, I do not duplicate the documentation and refer to the SciPy documentation 
	# http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html and 
	# http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html for 
	# further details.   
		  
	allFrag = fragSimTbl.keys()
	nfrag = len(allFrag)
	allFrag.sort()
	fragIdxTbl = {}

	fragDistArr = np.ones((nfrag,nfrag))
	for fi, frag in enumerate(allFrag):
		fragIdxTbl[frag] = fi
		fragDistArr[fi,fi] = 0.
	
	minDist = 1.
	maxDist = 0.  
	for frag1 in allFrag:
		fidx1 = fragIdxTbl[frag1]
		# self-distance = zero
		for frag2 in fragSimTbl[frag1]:
			fidx2 = fragIdxTbl[frag2]
			sim = fragSimTbl[frag1][frag2]
			# NB: convert to DISTANCE
			dist = 1. - sim
			if dist > maxDist:
				maxDist = dist
			if dist < minDist:
				minDist = dist
			fragDistArr[fidx1,fidx2] = dist
			# NB: need to redundantly put both symmetric distances or squareform() bitches
			fragDistArr[fidx2,fidx1] = dist

	distVec = distance.squareform(fragDistArr)
	
	try:
		clustering = fastcluster.linkage(distVec,method='ward')
	except Exception,e:
		print 'bldFragClusters: clustering exception',e
		import pdb; pdb.set_trace()
		
	# NB: sch.linkage doesn't support Ward method when "raw" observations not available
	# clustering = sch.linkage(distVec,method='average')
	
	flatClust = sch.fcluster(clustering, distThresh, 'distance')
	
	# clustRoot,clustTree = sch.to_tree(clustering, rd=True)
	
	# 2do:  better way to get most-central fragment for each cluster?
	# http://stackoverflow.com/questions/9362304/how-to-get-centroids-from-scipys-hierarchical-agglomerative-clustering
	
	clustTbl = defaultdict(list) # clustID -> [fragments]
	for fragIdx,cidx in enumerate(flatClust):
		clustTbl[cidx].append(allFrag[fragIdx])
	
	clustIndices = clustTbl.keys()
	
#	 cliqueSizes = [len(clustTbl[c]) for c in clustTbl]
#	 print 'bldFragClusters: NClique=%d %s' % (len(clustTbl),cliqueSizes)

	clustIndices.sort()
	centralFrag = {} # clusterIdx -> most central frag
	for cidx in clustIndices:
		ctrFrag = getCenter(clustTbl[cidx])
		centralFrag[cidx] = ctrFrag	

	fragClustTbl = {} # cliqueIdx -> (ctrFrag, [cliqueFrags] )				  
	for cidx in clustIndices:
		fragClustTbl[cidx] = (centralFrag[cidx], clustTbl[cidx])
			
	# print 'bldFragClusters: NFragment = %d minDist=%f maxDist=%f NClust=%d' % (nfrag,minDist,maxDist,len(clustIndices))

	return fragClustTbl

def loadLigPaths(inf):
	reader = csv.DictReader(open(inf))
	pathTbl = {} # ligand -> path
	for i,entry in enumerate(reader):
		# Lib,Ligand,Path
		pathTbl[ entry['Ligand'] ] = entry['Path']
	print 'loadLigPaths: %d ligand paths loaded' % (len(pathTbl))
	return pathTbl

def vdwpPP(patchList):
	'''convert list of ratoms, collapse shared chain,residue
	return as string: chain_residue:[intraResAtoms]
	'''
	
	resTbl = defaultdict(list)
	for vdwp in patchList:
		bits = vdwp.split('_') # chain,residue,ratom
		res = '+'.join(bits[:2]) 
		resTbl[res].append(bits[2])
	allRes = resTbl.keys()
	allRes.sort()
	outs = ''
	for ir,res in enumerate(allRes):
		if ir > 0:
			outs += ';'
		outs += '%s:[' % (res)
		raList = resTbl[res]
		raList.sort()
		for ir2, ra in enumerate(raList):
			if ir2 > 0:
				outs += ','
			outs += '%s' % ra
		outs += ']'
	
	return outs

def vdwPP2vdwpSet(vdwPP):
	'''return set of independent ratoms
	inverse for vdwpPP()
	'''
	
	# vdwPP ala 'B+029D:[H,N,CB];B+030D:[N,CG,CB,H,OD2];B+047I:[CB,CG2]'
	
	vdwRASet = set()
	resList = vdwPP.split(';')
	for res in resList:
		chainRes,ratomListStr = res.split(':')
		chain,resName = chainRes.split('+')
		ratomListStr = ratomListStr.replace('[','')
		ratomListStr = ratomListStr.replace(']','')
		ratomList = ratomListStr.split(',')
		for ratom in ratomList:
			vdwRASet.add( '_'.join([chain,resName,ratom]) )
	
	return vdwRASet

def vdwpRes(patchList):
	'''return list of receptor residues involved in vdwPatch
	'''
	resTbl = defaultdict(list)
	for vdwp in patchList:
		bits = vdwp.split('_')
		res = '_'.join(bits[:2])
		resTbl[res].append(bits[2])
	allRes = resTbl.keys()
	allRes.sort()
	
	return allRes

### top-level run commands
if __name__ == '__main__':

	HostName = socket.gethostname()
	if HostName == 'mgl0':
		print 'running on mgl0, good!'
		BaseDir = '/export/wcg/'
	
	elif HostName == 'mgl3':
		print 'running on mgl3, slow(:'
		BaseDir = '/mgl/storage/wcg/'
	
	elif HostName.startswith('hancock-vb'):
		print 'running on VBox-hancock'
		BaseDir = '/media/sf_sharedData/coevol-HIV/WCG/'	
	
	elif HostName.startswith('hancock'):
		print 'running local on hancock'
		BaseDir = '/Data/sharedData/coevol-HIV/WCG/'

	elif HostName.startswith('mjq'):
		print 'running local on mjq'
		BaseDir = '/home/Data/coevol-HIV/WCG/'

	else:
		print 
		sys.exit( ('unknown host %s' % (HostName)) )
	
	if len(sys.argv) < 2:
		sys.exit('FAAHA: missing runName argument?!')
	config.RunName = sys.argv[1]

	FocalExpts = None 
	LigPathTbl = None

	# NB: expt_160208 prefixed added!
	SummRptDir = BaseDir + 'anal/%s/'  % (config.RunName)

	ProcDir = BaseDir  + 'process2/'
	CrawlDir = BaseDir + 'crawl2/'


	############################################################
	
	print '<FAAHA config.RunName=%s begTime=%s>' % (config.RunName,datetime.datetime.now().strftime('%y%m%d_%H%M%S'))

	initFile = SummRptDir + 'config.ini'
	cfg = ConfigParser.ConfigParser()
	cfg.optionxform=str  # maintain options' case
	cfg.read(initFile)
	setConfigOptions(cfg, 'config')
	
	############################################################
	
	print '\thost=%s\n\tProcDir=%s\n\tCrawlDir=%s\n\tAnalDir=%s\n' % (HostName, ProcDir, CrawlDir, SummRptDir)

	## Top-level directories; created as needed below as part of each phase
	
	LowEDir = SummRptDir+'lowE/'
	InterTblDir = SummRptDir+'InterTable/'
	RLIFDir = SummRptDir+'RLIF_F/'

	FragSimDir = SummRptDir+ 'FragSim/'
	R2FDir = SummRptDir + 'R2F/'
	LigCoordDir = SummRptDir+'LigCoord/'
	R2FCDir = SummRptDir + 'R2FC/'
	PlotDir = SummRptDir + 'plots/'
	ArffDir = SummRptDir+'ARFF/'
	HIFDir = SummRptDir+'HIF/'
	VDWPDir = SummRptDir+'VDWPatch/'
	EliteDir = SummRptDir + 'EliteLig/'
	L2FDir = SummRptDir + 'L2F/'
	
	# NB: need to set PlotDir in config so plot.py() knows about it
	config.PlotDir = PlotDir

	# load mol2Idx
	config.Mol2Dir = None

# #	 ############################################################
		 
	print 'FAAHA: *PHASE 0: building exptTbl'
		  
	exptTblPkl0 = SummRptDir+'exptTblU_v0.pkl'
	if os.path.exists(exptTblPkl0):
		print '## Loading pre-pickled exptTbl_v0',exptTblPkl0
		exptTbl0 = cPickle.load(open(exptTblPkl0,'rb'))
	else:
		print '## <bldExptTbl>'
		exptTbl0 = bldExptTbl(SummRptDir+config.ExptFile)
		print '## </bldExptTbl>'
		print '## Saving pickled newExptTbl',exptTblPkl0
		cPickle.dump(exptTbl0, open(exptTblPkl0,'wb'))
	  
	if not os.path.isdir(LowEDir):
		print 'FAAHA: creating lowE directory',LowEDir
		os.makedirs(LowEDir)
	  
	exptTblPkl1 = SummRptDir+'exptTblU_v1.pkl'
	if os.path.exists(exptTblPkl1):
		print '## Loading pre-pickled exptTbl_v1',exptTblPkl1
		exptTbl = cPickle.load(open(exptTblPkl1,'rb'))

	else:
		print '## <FAAHA_expt>'
		exptTbl = FAAHA_expt(exptTbl0,FocalExpts,CrawlDir,SummRptDir)
		print '## </FAAHA_expt>'
		print '## Saving pickled exptTbl_v1',exptTblPkl1
		cPickle.dump(exptTbl, open(exptTblPkl1,'wb'))
	  
#	sys.exit('done building exptTbl')
	
#	 ############################################################
	# NB: subsequent phases all assume allExpt
		
	allExpt = exptTbl.keys()
	allExpt.sort()
			  
	############################################################

	print 'FAAHA: *PHASE 0-9: selectEliteSample'

	if not os.path.isdir(EliteDir):
		print 'FAAHA_0-9: creating EliteDir directory',EliteDir
		os.makedirs(EliteDir)

	for exptKey in allExpt:
		
		exptNo,prot,recept,site,lib = exptKey
		exptName = bldExptStr(exptKey)

		begTime = datetime.datetime.now()
		begTimeStr = begTime.strftime('%y%m%d_%H%M%S')
  
		print '<FAAHA_0-9 %s %s>' % (exptName,begTimeStr)

		exptData = exptTbl[exptKey]
		batchList = ranges2list( [(exptData['bstart'],exptData['bend'])] )
		qualRecept = exptData['receptor']

		eligf = EliteDir + exptName + '_eliteLig.csv'
		exptEliteDir = EliteDir + exptName + '/'
		
		# NB: presence of exptEliteLigFile used to flag prior compute
		# since exptEliteDir might have useful, previously getPDBQT.mgl_getADV_PDBQT() extracted PDBQTs from tarball
		if os.path.exists(eligf):
			print 'FAAHA_0-9: %s exptEliteLigFile exist; skipping' % (exptName)
			continue
		if not os.path.exists(exptEliteDir):
			os.makedirs(exptEliteDir)
		
		faahDir = CrawlDir + '%s/%s/' % (exptNo,qualRecept)

		uniqLigTbl, thresh, bestCandThresh = getThresh(exptKey,config.RunType,batchList,faahDir, \
													   config.Frac4Thresh,config.NCand,dcrit=config.Criterion)

		selectEliteDist(faahDir,exptKey,batchList,thresh,config.RabbleFrac,qualRecept)

		elapTime = datetime.datetime.now() - begTime
		print '</FAAHA_0-9 %s %s sec>' % (exptName,elapTime.total_seconds())

	############################################################
	 
	print 'FAAHA: *PHASE 1: building RLIF'
		   
	exptTblPkl1 = SummRptDir+'exptTblU_v1.pkl'
	exptTbl = cPickle.load(open(exptTblPkl1,'rb'))
		 
	if not os.path.isdir(InterTblDir):
		print 'FAAHA_1: creating InterTbl directory',InterTblDir
		os.makedirs(InterTblDir)
	   
	if config.vdwLigandAtom and not os.path.isdir(VDWPDir):
		print 'FAAHA_1: creating VDWPDir directory',VDWPDir
		os.makedirs(VDWPDir)
		   
	if not os.path.isdir(RLIFDir):
		# NB: presence of RLIFDir taken as indication this phase accomplished!
		print 'FAAHA_1: creating RLIF directory',RLIFDir
		os.makedirs(RLIFDir)

	for exptKey in allExpt:
		
		exptNo,prot,recept,site,lib = exptKey
		exptName = bldExptStr(exptKey)

		begTime = datetime.datetime.now()
		begTimeStr = begTime.strftime('%y%m%d_%H%M%S')
  
		print '<FAAHA_1 %s %s>' % (exptName,begTimeStr)

		loadNonZinc(exptName)
		
		if config.NBestLig == None:
			analRLIF(exptKey)
		else:
			analRLIF(exptKey,nbest=config.NBestLig)

		elapTime = datetime.datetime.now() - begTime
		print '</FAAHA_1 %s %s sec>' % (exptName,elapTime.total_seconds())
				   
#	 sys.exit('done building RLIF')

	############################################################
	  
	print 'FAAHA: *PHASE 1-1: compute HIF'
	if not os.path.isdir(HIFDir):
		print 'FAAHA_1-1: creating HIFDir directory',HIFDir
		os.makedirs(HIFDir)

	for exptKey in allExpt:
		
		exptNo,prot,recept,site,lib = exptKey
		exptName = bldExptStr(exptKey)
		
		hifFile = HIFDir+('%s.csv' % (exptName))
		if os.path.exists(hifFile):
			print 'FAAHA_1-1: %s HIF file exist; skipping' % (exptName)
			continue	
		
		begTime = datetime.datetime.now()
		begTimeStr = begTime.strftime('%y%m%d_%H%M%S')
  
		print '<FAAHA_1-1 %s %s>' % (exptName,begTimeStr)
		
		faahDir = CrawlDir + '%s/%s/' % (exptNo,recept)

		batchList = ranges2list( [(exptTbl[exptKey]['bstart'],exptTbl[exptKey]['bend'])] )

		frac4Thresh=config.Frac4Thresh			
		uniqLigTbl, thresh, bestCandThresh = getThresh(exptKey,config.RunType,batchList,faahDir,frac4Thresh, \
													   ncand=config.NCand,dcrit=config.Criterion)
		
		r2lf = InterTblDir + '%s_r2l.csv' % (exptName)
		r2lTbl = loadRLIF2LigTbl(r2lf)
				
		analBldHIFeatures(faahDir,exptName,batchList,thresh,r2lTbl,uniqLigTbl,hifFile)

		elapTime = datetime.datetime.now() - begTime
		print '</FAAHA_1-1 %s %s sec>' % (exptName,elapTime.total_seconds())

	# sys.exit('done building HIF')
			
	############################################################

	print 'FAAHA: *PHASE 2: compute lig2frag'
	if not os.path.isdir(L2FDir):
		print 'FAAHA_2: creating L2FDir directory',L2FDir
		os.makedirs(L2FDir)
	   
	for exptKey in allExpt:
		exptNo,prot,recept,site,lib = exptKey
		exptName = bldExptStr(exptKey)

		config.DropLigFile = SummRptDir + exptName + '_droppedLig.csv'
		dls = open(config.DropLigFile,'w')
		dls.write('ZINCID,Reason\n')
		dls.close()

		L2FPickleFile =   L2FDir + exptName + '_lig2Frag.pkl'
		if os.path.exists(L2FPickleFile):
			print 'FAAHA_2: %s L2F pkl exist; skipping' % (exptName)
			continue
					   
		begTime = datetime.datetime.now()
		begTimeStr = begTime.strftime('%y%m%d_%H%M%S')
  
		print '<FAAHA_2 %s %s>' % (exptName,begTimeStr)
				
		ligTbl, foo = provideEliteDist(exptName)
					   
		config.BindPPFile = L2FDir + exptName + '_lig2fragPP.txt'
			
		exptLigDir = EliteDir 
			
		if config.CheckMol2Lig:
			config.mol2Idx = bldDUDEMol2Idx(prot) # ligand -> mol2File
			
		l2fTbl = bldLig2frag(ligTbl,exptName,exptLigDir)
		   
		print 'L2FPickle dumping to',L2FPickleFile
		cPickle.dump(l2fTbl, open(L2FPickleFile,'wb'))
		
		fragDistTbl = fragDist(l2fTbl)
		fragDistFile = L2FDir + exptName + '_fragDist.csv'
		rptFragDist(fragDistTbl,fragDistFile)

		elapTime = datetime.datetime.now() - begTime
		print '</FAAHA_2 %s %s sec>' % (exptName,elapTime.total_seconds())
				
	# sys.exit( ('Done compute lig2frag' ) )
				  
	############################################################
  
	print 'FAAHA: *PHASE 3: build bldRLIF2frag'
	
	if not os.path.isdir(R2FDir):
		print 'FAAHA_3: creating RLIF2F directory',R2FDir
		os.makedirs(R2FDir)

	# nbest = 1000
	for exptKey in allExpt:  
		exptNo,prot,recept,site,lib = exptKey
		exptName = bldExptStr(exptKey)
  
		R2FPickleFile =   R2FDir + exptName + '_rlif2frag.pkl'
  
		if os.path.exists(R2FPickleFile):
			print 'FAAHA_3: %s R2F pkl exists; skipping' % (exptName)
			continue	
	   
		begTime = datetime.datetime.now()
		begTimeStr = begTime.strftime('%y%m%d_%H%M%S')
		print '<FAAHA_3  %s %s>' % (exptName,begTimeStr)
				
		ligTbl, foo = provideEliteDist(exptName)
	
		print 'FAAHA_3: Exp=%s: NLig=%d' % (exptName,len(ligTbl))
	   
		L2FPickleFile =   L2FDir + exptName + '_lig2Frag.pkl'
		if os.path.exists(L2FPickleFile):
			# print 'L2FPickle exists; using',L2FPickleFile
			l2fTbl = cPickle.load(open(L2FPickleFile,'rb'))
		else:
			print  'FAAHA: %s no L2FPickleFile found!? run bldLig2frag()' % (exptName)
			continue
							  
		print 'R2FPickleFile not found, building...',R2FPickleFile
		
		r2lf = InterTblDir + '%s_r2l.csv' % (exptName)
		l2rTbl = bldLig2RLIFTbl(r2lf)
		errf = R2FDir + exptName + '_r2f_err.txt'
		
		r2fTbl = bldRLIF2frag(exptName,ligTbl,l2fTbl,l2rTbl,errf,None) # RLIF -> frag -> [ (zincid,fragIdx,latomFull,fragInfo) ]
  
		print 'R2FPickleFile dumping to',R2FPickleFile
		cPickle.dump(r2fTbl, open(R2FPickleFile,'wb'))
		r2fFile = R2FDir + exptName + '_rlif2Frag.csv'
		rptR2FTbl(r2fTbl,r2fFile)
  
		elapTime = datetime.datetime.now() - begTime
		print '</FAAHA_3 %s %s sec>' % (exptName,elapTime.total_seconds())
	   
	# sys.exit( 'bldRLIF2frag done.' )
	 
		
#	 ############################################################
	 
	print 'FAAHA: *PHASE 5: bldR2FCtr'
		 
	if not os.path.isdir(R2FCDir):
		print 'FAAHA_5: creating R2FCDir directory',R2FCDir
		os.makedirs(R2FCDir)
			  
	for exptKey in allExpt:
	  
		exptName = bldExptStr(exptKey)
			 
		# NB: existence of BOTH fragClust, plot used to indicate phase already run
			 
		R2FCPickleFile = R2FCDir + exptName + '_r2fc.pkl'	
		if os.path.exists(R2FCPickleFile):
			print 'FAAHA Phase 5: fragClust already exist',exptName
			continue
   
		begTime = datetime.datetime.now()
		begTimeStr = begTime.strftime('%y%m%d_%H%M%S')
		print '<FAAHA_5  %s %s>' % (exptName,begTimeStr)
		 
		R2FPickleFile = R2FDir + exptName + '_rlif2frag.pkl'
		r2fTbl = cPickle.load(open(R2FPickleFile,'rb')) # RLIF -> frag -> [ (zincid,fragIdx,latomFull,fragInfo) ]
				  
		r2fcfile = R2FCDir + exptName + '_r2fc.csv'
		
#		 hifFile = HIFDir+('%s.csv' % (exptName))
#		 rlifInfoTbl = loadRLIFInfo(hifFile)
		
		r2fcTbl = bldR2FCtr(exptName,r2fTbl,None,r2fcfile)
		   
		cPickle.dump(r2fcTbl, open(R2FCPickleFile,'wb'))
   
		elapTime = datetime.datetime.now() - begTime
		print '</FAAHA_5 %s %s sec>' % (exptName,elapTime.total_seconds())
	 
	# sys.exit( 'bldRLIF2fragCtr  done.' )
	
	if not os.path.isdir(ArffDir):
		print 'FAAHA_6: creating ArffDir directory',ArffDir
		os.makedirs(ArffDir)
			
	print 'FAAHA: *PHASE 6: build ARFF'
   
	# NB: need to set PlotDir in config so plot.py() knows about it
	if not os.path.isdir(config.PlotDir):
		print 'FAAHA_6: creating PlotDir directory', config.PlotDir
		os.makedirs(config.PlotDir)
	 
	for exptKey in allExpt:
		exptNo,prot,recept,site,lib = exptKey
		exptName = bldExptStr(exptKey)
  
		begTime = datetime.datetime.now()
		begTimeStr = begTime.strftime('%y%m%d_%H%M%S')
  
		print '<FAAHA_6 %s %s>' % (exptName,begTimeStr)
		
		arrfFile = ArffDir+('%s_' % (exptName))
		if not config.FragQualAttrib:
			arrfFile += 'noFrag_'
		if not config.VDWPAttrib:
			arrfFile += 'noVDWP_'
		arrfFile += 'sp.arff' 
		
		print 'FAAHA_6: ARFF file:',arrfFile
		
		if os.path.exists(arrfFile):
			print 'FAAHA_6: arrf already exist',exptName
			continue
						   
		ligTbl, foo = provideEliteDist(exptName)
		activeIdxSet = provideEliteTrain(exptName)
			
		cummTrueMaxE, fracNeg = analTrueEnergy(exptName, ligTbl, activeIdxSet)
		print 'FAAHA_6: %s cummTrueMaxE=%f FracNeg=%f (%d)' % \
			(exptName,cummTrueMaxE,fracNeg,fracNeg*len(ligTbl))
			
		R2FCPickleFile = R2FCDir + exptName + '_r2fc.pkl'
		r2fcTbl = cPickle.load(open(R2FCPickleFile,'rb'))
		
#		 lcFile = LigCoordDir + '%s_ligCoord.pkl' % (exptName)
#		 print 'Loading ligCoordPickle...',lcFile
#		 ligCoordTbl = cPickle.load(open(lcFile,'rb'))

		R2FPickleFile = R2FDir + exptName + '_rlif2frag.pkl'
		r2fTbl = cPickle.load(open(R2FPickleFile,'rb')) # RLIF -> frag -> [ (ligIdx,fragIdx,latomFull,fragInfo) ]

		lig2SpArff(exptName,ligTbl,r2fcTbl,r2fTbl,activeIdxSet,arrfFile,\
				ethresh=None,useFragQual=config.FragQualAttrib,useVDWP=config.VDWPAttrib,splitRAtomAttrib=config.SplitRAtomAttrib)

		elapTime = datetime.datetime.now() - begTime
		print '</FAAHA_6 %s %s sec>' % (exptName,elapTime.total_seconds())
			
	# sys.exit('done creating ARRF')
  
