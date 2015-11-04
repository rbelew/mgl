''' crawl_ADV: utilities to crawl ADVina FAAH files
    was part of get_ADInfo

@version 1.0
@date on 30 Aug 14
@author: rbelew@ucsd.edu, dsantiago@scripps.edu
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

# from AutoDockTools.HelperFunctionsN3P import pathToList, getLines, percent
from string import strip

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
#                'GLN':'Q', 'ASN':'N', 'SER':'S', 'ASX':'B', 'GLX':'Z',
#                'PHE':'F', 'TRP':'W', 'TYR':'Y',
#                'GLY':'G', 'ALA':'A', 'ILE':'I', 'LEU':'L', 'CYS':'C',
#                'MET':'M', 'THR':'T', 'VAL':'V', 'PRO':'P' }

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
InterRE = r'(.*):.+():([A-Z]+[0-9]*)~~(.*):(.+):(.+)'

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
            if itype=='vdw':
                for inter in ligData[itype]:
                    # d:<0>:O3~~B:ARG57:N
                    # :LEU63:CD1 -->actual example, matches InterVDWREPat \dns
                    m = InterVDWREPat.match(inter)
                    try:
                        (rchain,raa,ratom) = m.groups()
                        rinterDict[itype].append( (rchain,raa,ratom) )
                    except:
                        print 'reducePlusInterDict: bad vdw string?!',inter
            elif itype=='ppi' or itype=='tpi':
                for inter in ligData[itype]:
                    # d:<0>:O3~~B:ARG57:N
                    m = InterPiREPat.match(inter)
                    try:
                        (rchain,raa,rcenter,ligcenter) = m.groups()
                        rinterDict[itype].append( (rchain,raa,rcenter,ligcenter) )
                    except:
                        print 'reducePlusInterDict: bad ppi/tpi string?!',inter
            else:
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
    if lines[0].startswith("USER    ADVS_Vina_result>"):
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
    for l in lines:
        if l.startswith(srcPrefix):
            src = l[len(srcPrefix):].strip()
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
    
def parseADPDBQT_ADV(f):
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
                   
    return ligData
        

# simplified, here, means that this should work for all batches after ??? Exp.# ???
# because naming was simplified (by dns)
# ADVsimple = r'fahv.x([a-zA-Z0-9]+)(_[A-Z0-9]*)*_(ZINC[0-9]+)(_[0-9]*)*_([0-9]+)_out_Vina_VS.pdbqt'
ADVsimple = r'fahv.x([a-zA-Z0-9]+)_(.+)_([0-9]+)_out_Vina_VS.pdbqt'
ADVsimplePat = re.compile(ADVsimple, re.IGNORECASE)

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
    #    fahv.x3kf0A_ZINC01569654_1113915765_out_Vina_VS.pdbqt

    procList = glob.glob(tmpDir+'/FAHV*/fahv.*_out_Vina_VS.pdbqt')
    if verbose:
        print 'visit_ADV_tgz: NTGZ=',len(procList)
        
    for isd,procPath in enumerate(procList):
        # fahv.x3kf0A_ZINC00145439_2057149382_out_Vina_VS.pdbqt
        procBits = os.path.split(procPath)

        # NB: don't need to capture batchNo; visitRpt_ADV_tgz() has it!
        # /tmp/tmpmmBQ7D/FAHV_x3kf0A_0124412_processed
        # dirBits = procBits[0].split('_')
        # batchNo = int(dirBits[2])

        procf = procBits[1]
        
        # lbpos = procf.find('_')
        # rbpos = procf.rfind('_')
        # assert (lbpos != -1 and rbpos != -1 and rbpos > lbpos), 'visit_ADV_tgz: bad procf?! %s' % (procf)
        # ligand = procf[lbpos+1:rbpos]
        
        mpath = ADVsimplePat.match(procf)
        (receptor,ligand,workNo) = mpath.groups()
        
        # import pdb; pdb.set_trace()
        
        ###-------
        ligData = parseADPDBQT_ADV(procPath)
        ###-------
        
        if not(ligData):
            print 'visit_ADV_tgz: invalid ADV file?!',procf, tgzPath
            continue

        dk = (exptname,receptor,ligand)
        
        if dk in dataTbl:
            print 'visit_ADV_tgz: dup dataKey?!',dk
            continue

        dataTbl[dk] = ligData

                                
    shutil.rmtree(tmpDir)
    # print 'visit_ADV_tgz: done.',len(dataTbl)
    
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


def visitRpt_ADV_tgz(tgzPath,recon,batchTbl,outdir,tocs,exptname,batchNo,verbose):
    ''' get info, place in this table
    	info is extracted from enhanced pdbqt files inside tarball
    '''	

    dataTbl = visit_ADV_tgz(tgzPath,exptname,recon,verbose)
    if recon:
        print 'visitRpt_ADV_tgz: Recon-only; no reporting'
        return len(dataTbl)
        
    summf  = outdir+exptname+'/summ/ADV_summ_%07d.csv' % (batchNo)
    interf = outdir+exptname+'/inter/ADV_inter_%07d.json' % (batchNo)
    rptData_ADV(dataTbl,summf,interf,exptname,batchNo)
    
    #tocStr = '%s,%05d,%d,%s' % (exptname,ntgz,len(dataTbl),tgzPath)
    tocStr = '%s,%d,%d,%s' % (exptname,batchNo,len(dataTbl),tgzPath)
    if verbose:
        print 'visitRpt_ADV_tgz toc:',tocStr
            
    tocs.write(tocStr+'\n')
    # fun2watch! toc
    tocs.flush(); os.fsync(tocs.fileno())
    
    return len(dataTbl)

def mglTop_visit_ADV(ADV_topDir,outdir,exptList=None,recon=False,verbose=False):
    'recon stops after opening, parsing one file in first tgz'
    
    if exptList:
        print 'mglTop_visit_ADV: Explicit experiment list %s' % (str(exptList))
    else:
        crawlPat = ADV_topDir+'/Exp*'
        print 'mglTop_visit_ADV: Full crawl of %s' % (crawlPat)
        exptList = [os.path.split(exptPath)[1] for exptPath in glob.glob(crawlPat) ]
        
    print 'mglTop_visit_ADV: NExperiments=',len(exptList)
    
    if verbose:
        print 'mglTop_visit_ADV: **Verbose output'

    if recon:
        print 'mglTop_visit_ADV: **Reconnaissance sniffing only!'

    totParse = 0
    for ie, exptname in enumerate(exptList):

        startTime = datetime.datetime.now()
        print 'mglTop_visit_ADV: %s starting %s' % (exptname,startTime.strftime("%y%m%d %H:%M:%S"))

        exptPath = ADV_topDir+exptname
        outPath = outdir+exptname
        
        if verbose:
            print ' *',ie,exptPath
            
        try:
            tst = open(outPath+'/tst.csv','w')
        except:
            print 'mglTop_visit_ADV: creating ExptOutput directory', (outPath)
            os.makedirs(outPath)
        tocf = outPath+'/ADV_toc.csv'
        tocs = open(tocf,'w')
        #tocs.write('NTGZ,Data,Path\n')
        tocs.write('Experiment, Batch, Data, Path\n')

        try:
            tst = open(outPath+'/summ/tst.csv','w')
        except:
            print 'mglTop_visit_ADV: creating summ directory',outPath+'/summ'
            os.makedirs(outPath+'/summ')
        try:
            tst = open(outPath+'/inter/tst.csv','w')
        except:
            print 'mglTop_visit_ADV: creating inter directory',outPath+'/inter'
            os.makedirs(outPath+'/inter')

        batchTbl = {}

        exptSubList = glob.glob(exptPath+'/Results_*')   # Results_* level dirs containing
#                                                     FAHV_<rec>_<batch_num>_processed.tgz files
        print 'mglTop_visit_ADV: NSubExperiments=',len(exptSubList)
        for ise, exptSubPath in enumerate(exptSubList):
            if verbose:
                print ' **',ie,ise,exptSubPath

            tgzList = glob.glob(exptSubPath+'/*.tgz') 
            print 'mglTop_visit_ADV: NTGZ=',len(tgzList)
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
                    print 'Attempting to analyze',tgznow 
                                    
                # process all enhanced pdbqt files in the processed batch file: ...
                #
                nparse = visitRpt_ADV_tgz(tgzPath,recon,batchTbl,outdir,tocs,exptname,batchNo,verbose)
                #
                totParse += nparse
                
        endTime = datetime.datetime.now()
        elapTime = endTime-startTime
        print 'mglTop_visit_ADV: %s done. TotParse=%d NSec=%d' % (exptname,totParse,elapTime.seconds)
            
        tocs.close() # for each experiment directory
        
    print 'mglTop_visit_ADV: TotParse=',totParse

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
#                                  ...Results_x1FB7-1w0n0-AAsh25/PDB017_x1FB7-1w0n0-AAsh25_vina-out_Vina_VS.pdbqt
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

import subprocess
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
    

def Local_Top_visit_ADV(fldir,outdir,exptname,sharedRecept,fileList,verbose=False):
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
        
        if sharedRecept == None:
            if RunName == 'SAMPL4':
                pathName,basename = os.path.split(path)
                pbits2 = pathName.split('/')
                receptName = pbits2[-2]
#                 print receptName
#                 receptorPath = pathName + '/' + receptName + '.pdbqt'
            else:
                print 'Local_Top_visit_ADV: missing receptor?!',RunName
                receptName = '??'
                continue
        else:
            fname = pathBits[-1]
            ppos = fname.find('.')
            receptName = fname[:ppos]
#             receptPath = sharedRecept

        dk = (exptname,receptName,ligand)

        # assert dk not in allPaths, 'Local_Top_visit_ADV: dup ligand?! %s %s' % (lib,ligand)
        if dk in allPaths:
            print 'Local_Top_visit_ADV: dup dataKey?! %s' % (dk,)
            ndup += 1
            
        allPaths[dk] = path
        relPath = path.replace(ProcDir,'')
        outs.write('%s,%s,%s,%s,%s\n' % (exptname,lib,receptName,ligand,relPath))
                
        if dk in dataTbl:
            print 'Local_Top_visit_ADV: dup dataKey?!',dk
            continue

        dataTbl[dk] = ligData
                               
        totParse += 1
    
    outs.close()        
    print 'Local_Top_visit_ADV: TotParse=%d NDup=%d' % (totParse,ndup)
    
    summf  = outdir+ ('ADV_summ_%07d.csv' % (1))
    interf = outdir+ ('ADV_inter_%07d.json' % (1))
    rptLocalData_ADV(dataTbl,summf,interf)
    
## mgl3 crawl of AD                                                                                                                                          
# import socket
# if socket.gethostname() == 'mgl0':
#     print 'running on mgl0, good!'
#     ADV_topDir = '/export/wcg/processed/'
#     outdir = '/export/wcg/crawl/test/'
#
# elif socket.gethostname() == 'mgl3':
#     print 'running on mgl3, slow(:'
#     ADV_topDir = '/mgl/storage/wcg/processed/'
#     outdir = '/mgl/storage/wcg/crawl/test/'
#
# AD_topDir = '/Data/sharedData/coevol-HIV/WCG/subsets/FAHV_tst_140829_2/'
# outdir = '/Data/sharedData/coevol-HIV/WCG/summRpts/tst/FAHV_140829/'
# exptList = ['Exp81','Exp82', 'Exp84'] # CB, AS, MB; smallest
#
# mglTop_visit_ADV(ADV_topDir, outdir, exptList)

# arg string ala:
# 140905
# /mgl/storage/wcg/processed/ /mgl/storage/wcg/crawl/test/  --exptList "['Exp72']" --verbose
# 140906
# /mgl/storage/wcg/processed/ /mgl/storage/wcg/crawl/test/  --exptList "# arg string ala:
# 140905
# /export/wcg/processed/ /export/wcg/crawl/test/  --exptList "['Exp81','Exp82','Exp84']" --verbose

# 141030: local rikPR253
# /Data/sharedData/coevol-HIV/WCG/processed/rikPR253/processed /Data/sharedData/coevol-HIV/WCG/crawl/rikPR253 --verbose

parser = argparse.ArgumentParser(description='crawl_ADV arguments')
parser.add_argument('ADV_topDir',type=str,help='Path to crawling rootdir')
parser.add_argument('outDir',type=str,help='Path to directory to contain result files')
parser.add_argument('--exptListStr', action="store",help='list of subset exptDir to crawl(string)')
parser.add_argument("--verbose",  action="store_true",help="increase output verbosity")
parser.add_argument("--recon",  action="store_true",help="Reconnaissance sniffing only")

# if __name__ == '__main__': 
#     
#     args, unknown = parser.parse_known_args()
#     if args.verbose:
#         print 'crawl_ADV: arguments'
#         # NB: args is a Namespace object; 
#         argsDict = vars(args)
#         for k,v in argsDict.items():
#             print '\t%s = %s' % (k,v)
#     
#     if len(unknown)>0:
#         print 'crawl_ADV: huh?! Unkown arguments=', unknown
#         assert False # can't do break or return here!
#     
#     if args.exptListStr:
#         exptList = eval(args.exptListStr)
#     else:
#         exptList = None
#         
#     # print '# PROFILING run!!'
#     # import cProfile
#     # cProfile.run(('mglTop_visit_ADV("%s","%s",%s,verbose=%s)' % (args.ADV_topDir, args.outDir, exptList, args.verbose)), \
#     #              args.outDir+'ADV_mgl0_profile.txt')
# 
#     # mglTop_visit_ADV(args.ADV_topDir, args.outDir, exptList, verbose=args.verbose)
# 
#     Local_PR253_Top_visit_ADV(args.ADV_topDir, args.outDir, verbose=args.verbose)

### top-level run commands
if __name__ == '__main__':

    RunName = 'In-LEDGF' # 'focusedLib'

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

    else:
        print 
        sys.exit( ('unknown host %s' % (HostName)) )

    ProcDir =    BaseDir + 'processed/%s/'  % (RunName)
    CrawlDir = BaseDir + 'crawl/%s/'  % (RunName)
    SummRptDir = BaseDir + 'anal/%s/'  % (RunName)

    if RunName == 'SAMPL4':
        ProcDir += 'LEDGF/'
        dockfsuffix = '*_VS.pdbqt'
        cmdList = ['find', ProcDir, '-name', dockfsuffix]
        outstr = subprocess.check_output(cmdList)
        fileList = outstr.split()
        print 'crawl_ADV: %s NFiles=%d' % (RunName,len(fileList))
        receptor = None
        
    elif RunName == 'focusedLib':
        dockfsuffix1 = '*_VS.pdbqt'
        dockfsuffix2 = '*.VS.pdbqt'

        cmdList1 = ['find', ProcDir, '-name', dockfsuffix1]
        outstr1 = subprocess.check_output(cmdList1)
        fileList1 = outstr1.split()
    
        cmdList2 = ['find', ProcDir, '-name', dockfsuffix2]
        outstr2 = subprocess.check_output(cmdList2)
        fileList2 = outstr2.split()
    
        fileList = fileList1 + fileList2
        
        print 'crawl_ADV: %s NFiles=%d (_VS=%d .VS=%d)' % (RunName,len(fileList),len(fileList1),len(fileList2))

        receptor = '/Data/sharedData/coevol-HIV/WCG/processed/focusedLib/lib_v1.1/CCDKF115_dimer.pdbqt'

    if not os.path.isdir(CrawlDir):
        print 'crawl_ADV: creating CrawlDir directory',CrawlDir
        os.makedirs(CrawlDir)

    Local_Top_visit_ADV(ProcDir,CrawlDir,RunName,receptor,fileList)
