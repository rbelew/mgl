''' analFAAH: top-level full-workflow module
Created on Jun 23, 2015

@author: rik
'''

 
from collections import defaultdict
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
import sys
  
import numpy as np

import scipy.cluster.hierarchy as sch
# NB, using scipy's spatial distance (not SKLearn's metrics)
import scipy.spatial.distance as distance
from scipy.sparse.dok import dok_matrix
 
import matplotlib
matplotlib.use('Agg') # ASSUME no windowing
# Although many examples use pylab, it is no longer recommended.                                                         
# import pylab as p
import matplotlib.pyplot as p
 
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
 
import pandas as pd

## other analFAAH modules
import config
import crawl_ADV
import plot

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

# obsolete?
LL2HIFDir = ''
SAMPLHIFDir = ''
StatsDir = ''
Best2RLIFDir= ''
ActiveDir = ''

InterTypes = ('hba', 'hbd', 'mtl','ppi','tpi','vdw')
BinaryITypes = ('hba', 'hbd', 'mtl','ppi','tpi')

FE_coeff_tors = 0.2983 # cf AD4.1_bound.dat

# 1 Oct 13: HACK: repair parse of experiment names including ligand suffix (:

OddExptRE = '^([0-9]+_)(.*)'
OddExptREPat = re.compile(OddExptRE)

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


## utilities
def basicStats(l):
    "Returns avg and stdev"
    if len(l) == 0:
        return(0.,0.)

    sum = 0
    for n in l:
        sum += n
    avg = float(sum) / len(l)

    sumDiffSq = 0.
    for n in l:
        sumDiffSq += (n-avg)*(n-avg)

    stdev = math.sqrt(sumDiffSq) / float(len(l))
    return (avg,stdev)

def clean_axis(ax):
    """Remove ticks, tick labels, and frame from axis"""
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

def frange(x, y, jump):
    while x < y:
        yield x
        x += jump

def entropy(dist):
    'compute (log2 ==> bits) entropy over distribution'

    h = 0.
    tot = sum(dist)
    if tot == 0:
        return h
    
    for n in dist:
        if n > 0:
            nf = float(n)
            p = nf/tot
            h -= p * math.log(p,2)
    return h

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

def bldRanges(l):
    "Identify contiguous ranges of integers; tuples iff more than one element in range"
    ll = l[:]  # don't want to apply in-place sort to list passed in!
    ll.sort()
    rangeList = []
    prev = 0
    begr = 0
    for bk in ll:
        if bk != prev+1:
            if prev != 0:
                if begr==prev:
                    rangeList.append( begr )
                else:
                    rangeList.append( (begr,prev) )
            begr = bk
        prev = bk
    if begr==prev:
        rangeList.append( begr )
    else:
        rangeList.append( (begr,prev) )
    return rangeList

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

def findex(frange,val):
    eps = 1e-6
    for i,v in enumerate(frange):
        if abs(frange[i]-val) < eps:
            return i 
    return -1
        

def avgTbl(tbl):
    'return average value over tables entries'

    avg =  float( sum(  [ tbl[e] for e  in tbl.keys() ] )) / len(tbl)
    return avg

def touchTimeDiff(f1,f2):
    '''given two files (earlyFile, lateFile), return seconds between their touch times
    '''
    t1 = os.stat(f1).st_mtime
    t2 = os.stat(f2).st_mtime
    
    return t2 - t1

def collectTimes(exptTbl,outf):
    # compare across different experiments
    interExpt = [   ['lowE', LowEDir, '_lowE.csv'],
                    ['interact', InterTblDir, '_interTbl.pkl'],
                    ['RLIF', RLIFDir, '_rlif.csv'],
                    ['l2f', L2FDir, '_lig2Frag.pkl'],
                    # ['fragSim', FragSimDir, '_fragSim.csv'],
                    ['R2F', R2FDir, '_rlif2frag.pkl'],
                    ['ligCoord', LigCoordDir, '_ligCoord.pkl'],
                    ['R2FC', R2FCDir, '_r2fc.pkl']
                    # ['arrf', ArffDir, '.arff']   # 2do: allow optional files
                    ]

    # compare different files within same experiment
    intraExptPairs = [] # [ ['fragClust', [FragSimDir, '_fragSim.csv'], [FragSimDir, '_fragClust.csv'] ]  ]

    allLbls = ['lowE', 'interact', 'RLIF', 'l2f','R2F', 'ligCoord', 'R2FC'] # 'arrf'

    allExpt = exptTbl.keys()
    allExpt.sort()
    
    exptNameList = [bldExptStr(expt) for expt in allExpt]    
   
    timeTbl = {ename: {} for ename in exptNameList}
    
    for ie,expt in enumerate(allExpt):
        # exptNo,prot,recept,site,lib = expt
        exptName = bldExptStr(expt)
        
        for traPair in intraExptPairs:
            lbl, fspec1, fspec2 = traPair
            f1 = fspec1[0] + exptName + fspec1[1]
            f2 = fspec2[0] + exptName + fspec2[1]
            timeTbl[exptName][lbl] = touchTimeDiff(f1, f2)
            
        # no previous-expt for comparision with first experiment(:
        if ie==0:
            continue
        
        for terspec in interExpt:
            lbl, fdir, fsuff = terspec
            f1 = fdir + exptNameList[ie-1] + fsuff
            f2 = fdir + exptName + fsuff
            timeTbl[exptName][lbl] = touchTimeDiff(f1, f2)

    outs = open(outf,'w')
    outs.write('Expt,'+ ','.join(allLbls) + '\n')
    for exptName in exptNameList:
        outs.write('%s' % (exptName))
        for lbl in allLbls:
            if lbl in timeTbl[exptName]:
                outs.write(',%d' % (timeTbl[exptName][lbl]))
            else:
                outs.write(', ')
        outs.write('\n')
    outs.close()

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
        
        if entry['experiment'].startswith('# '):
            continue
        
        # HACK:      
        # replace any '_' in receptor, site, library with '-'
        # to avoid interactions with underbar concatenations later
        
        exptData = {}

        try:
            exptNo = entry['experiment']
            prot = entry['protein'].replace('_','-')
            recept = entry['receptor'].replace('_','-')
            recept = entry['receptor'].replace('.pdbqt','')
            site = entry['site'].replace('_','-')
            lib = entry['library'].replace('_','-')
            
            exptData['bstart']  = int(entry['bstart'])
            exptData['bend']    = int(entry['bend'])
            exptData['sys']     = entry['program']
        except Exception,e:
            print 'bldExptTbl: bad entry?! %s %s' % (e,entry)
            continue
        
        exptData['lib']     = lib
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
    
def analFAAH_expt(exptTbl,exptList,faahDirParent,outDir,frac4Thresh=0.02,ncand=1000,dcrit='energy'):
    ''' Identify top ncand and top frac4Thresh for each experiment
    
        first pass sorts all energies, identifies threshold, 
        outputs best ncand  candidates ACROSS ALL EXPERIMENTS to best_*.csv
        two thresholds computed: threshold=thresh@frac4Thresh; bestCandThresh=thresh@ncand
        NB: exptTbl augmented with threshold!
        built from analFAAH, but for only focal experiments
    '''
                
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
            
    print 'analFAAH_expt: exptList=%s NExpt=%d' % (exptList,len(expt2do))

    candf = outDir+('bestLig_%d.csv' % ncand)

    cands = open(candf,'w')
    cands.write('Expt,Batch,Ligand,E,FullExpt\n')
    
    newExptTbl = {}
    expt2do.sort()
    for exptKey in expt2do:
        exptNo,prot,exptRecept,site,lib = exptKey
        exptData = exptTbl[exptKey]
        exptName = bldExptStr(exptKey)
        if config.RunName.startswith('SAMPL4') or \
            config.RunName.startswith('focusedLib'):
            batchList = [1]
        else:
            batchList = ranges2list( [(exptData['bstart'],exptData['bend'])] )
            
        sys = exptData['sys']
        lib = exptTbl[exptKey]['lib']

        nbatch = len(batchList)
        
        if config.RunName.startswith('SAMPL4') or \
            config.RunName.startswith('focusedLib'):
            faahDir = faahDirParent
        else:
            faahDir = faahDirParent + 'Exp%s/' % (exptNo)
            
        thresh, bestCandThresh = getThresh(exptKey,sys,batchList,faahDir,frac4Thresh,ncand,dcrit)
                
        newExptTbl[exptKey] = exptData.copy()
        
        # NB: exptTbl augmented with threshold!
        newExptTbl[exptKey]['thresh'] = thresh
        # 150930  keep ncand thresh info too
        newExptTbl[exptKey]['ncand'] = ncand
        newExptTbl[exptKey]['bestCandThresh'] = bestCandThresh

        ## 2d pass

        ## NonZincLigTbl is experiment-specific
        ## needs to be kept for reloading along with pickeled
        config.NonZincLigTbl = {} # ligand name -> nonzincID 
        config.NZIdx2LigTbl = {} # nonzincID -> ligand name
        config.NNonZincLig = 0

        outs = open(LowEDir+('%s_lowE.csv' % exptName),'w')
        # cf. mglTop_visit_AD() for header line
        fldNameLine = 'Batch,Ligand,E,Eff,DCVal,Nvdw,Ninter\n'
        outs.write(fldNameLine)
    
        nrcd = 0
        nlowE = 0
        ndup = 0
        lowETbl = {}
        ligTbl = {}
        
        for isf,batchNo in enumerate(batchList):
            # 140604
            if config.RunName.startswith('SAMPL4') or \
                config.RunName.startswith('focusedLib'):
                summPath = faahDir + 'ADV_summ_0000001.csv' 
            else:
                summPath = faahDir+('summ/%s_summ_%07d.csv' % (sys,batchNo))
                
            try:
                inStr = open(summPath)
            except Exception,e:
                # ASSUME counted as nmissf above
                # print "analFAAH_expt: can't open2 %s?!" % (summPath)
                continue
         
            for il,line in enumerate(inStr.readlines()):
                if il == 0:
                    # Expt,Recept,Ligand,E,Eff,Nvdw,Ninter
                    # written by get_ADInfo.rptData()
                    continue
                nrcd += 1
                try:
                    flds =line[:-1].split(',')

                    if config.RunName.startswith('SAMPL4') or \
                        config.RunName.startswith('focusedLib'):
                        (ligRecept,ligand,e,eff,nvdv,ninter) = flds
                    else:
                        # ['Exp32', '13651', 'xEyeSiteXtl5NI', 'ZINC00039702', '-4.63', '-0.421', '28', '4']
                        (expt,batch,ligRecept,ligand,e,eff,nvdv,ninter) = flds

                    if ligRecept != exptRecept:
                        continue
                        
                    if dcrit == 'energy':
                        dval = float(e)
                    elif dcrit == 'ligEff':
                        dval = float(eff)
                    elif dcrit == 'ligEnth':
                        # Assume LigDOFTbl already loaded
                        if ligand not in LigDOFTbl:
                            print 'analFAAH_expt: lig missing from LigDOFTbl?!',ligand
                            continue
                        dof = LigDOFTbl[ligand]
                        natom = float(e) / float(eff)
                        natomi = int(natom)
                        assert natomi-natom > 1e-2, 'huh?'
                            
                        dval = (float(e) + (dof * FE_coeff_tors)) / natom

                except Exception, e:
                    print "analFAAH_expt: bad line(2) %s %d %s?!" % (summPath,il,e)
                    break # don't try to read other lines from inStr
                
                # Pass2: need to make this thresh criterion-sensitive, too
                
                if (dcrit == 'energy' and abs(dval) > TooBigE) or \
                    (dcrit == 'ligEff' and abs(dval) > TooBigLigEff) or \
                    (dcrit == 'ligEnth' and abs(dval) > TooBigLigEnth):
                    continue

                ligIdx = normLigand(ligand)
                
                if config.RunName.startswith('DUDE'):
                    ligand2 = ligIdx2zinc(ligIdx)
                else:
                    ligand2 = ligIdx2zinc(ligIdx)
                if ligand2 != ligand:
                    # NB: two ligand cleanups!
                    if not (ligand.endswith('.VS') or ligand.find('pras') != -1):
                        print 'analFAAH_expt: bad ligand indexing?!', ligand, ligIdx, ligand2
                                                         
                ## dupes:  just count them
                if ligIdx in ligTbl:
                    ndup += 1
                else:
                    ligTbl[ligIdx] = 1
    
                ## thresh
                
                if dval < thresh:
                    nlowE += 1
                    ## NB: all 3 vals (e,eff,dval) written to lowE file; why not!
                    newline = '%d,%s,%s,%s,%s,%s,%s\n' % (batchNo,ligand,e,eff,dval,nvdv,ninter)      
                    outs.write(newline)
                if dval < bestCandThresh:
                    cands.write('%s,%d,%s,%f,%s\n' % (exptNo,batchNo,ligand,dval,exptName))
                
            inStr.close() # eo-batch file
            
        nlig = len(ligTbl)
            
        outs.close()
        
        if config.NNonZincLig>0:
            # cf. bldNonZincIdx()
            if config.RunName.startswith('DUDE'):
                nzligFile = LowEDir + '%s_nonZincLig.csv' % (exptName)
            else:
                nzligFile = LowEDir + 'nonZincLig.csv'
                
            print 'analFAAH_expt: Writing %d NonZinc ligands to %s' % (config.NNonZincLig,nzligFile)
            allLig = config.NonZincLigTbl.keys()
            allLig.sort()
            outs = open(nzligFile,'w')
            outs.write('Ligand,NZIdx\n')
            for lig in allLig:
                outs.write('%s,%d\n' % (lig,config.NonZincLigTbl[lig]))
            outs.close()

        try:
            avgRcd = float(nrcd)/nbatch
        except:
            avgRcd = 0
            
        try:
            fracLowE = float(nlowE)/nrcd
        except:
            fracLowE = 0
#        print 'analFAAH: %s NBatch=%d NRcd=%d NLigand=%d NLowE=%d NDup=%d AvgRcd/File=%f FracLowE=%f'  % \
        print '%s,%s,%d,%d,%d,%d,%d,%f,%f'  % \
                (exptName,lib,nbatch,nrcd,nlig,nlowE,ndup,avgRcd, fracLowE)
        
        # eo-expt
    cands.close()
    return newExptTbl

def getEvalLigands(inf,keyfld='PrTrue'): 
    'create ligand list from ligPredict.csv sorted by descending PrTrue'
    
    reader = csv.DictReader(open(inf))        
    ligTbl = {}
    for i,entry in enumerate(reader):
        # Ligand,Actual,Predict,PrTrue,Err,E,FPRate,TPRate
        ligTbl[entry['Ligand']] = float(entry[keyfld])
        
    ligands = ligTbl.keys()
    ligands.sort(key=lambda k: ligTbl[k], reverse= True)
    return ligands

def getClassifLigands(inf,keyfld='PrTrue'): 
    '''create (ligand,active,classLbl,errp) list from ligPredict.csv sorted by descending PrTrue
    for use by evalClassif.getCommonClassifLig()
    '''
    
    reader = csv.DictReader(open(inf))        
    ligTbl = {}
    for i,entry in enumerate(reader):
        # Ligand,Actual,Predict,PrTrue,Err,E,FPRate,TPRate
        ligTbl[entry['Ligand']] = (float(entry[keyfld]),entry['Actual'],entry['Predict'],entry['Err'])
        
    ligands = ligTbl.keys()
    ligands.sort(key=lambda k: ligTbl[k][0], reverse= True)
    
    ligList = []
    for lig in ligands:
        ligList.append( (lig,ligTbl[lig][1],ligTbl[lig][2],ligTbl[lig][3]) )
    return ligList

def getThresh(exptKey,sys,batchList,faahDir,frac4Thresh,ncand,dcrit):

    ## first-pass: set threshold

    ## 141219: selection criterion made a parameter
    # this requires mods to BOTH passes!
    # NB: dval used as generic value for Pass1, then retain allE(), variable names
    # "E" name is no longer correct, but doesn't hurt anything!

    exptNo,prot,exptRecept,site,lib = exptKey

    allE = []
    noddE = 0
    nmissf = 0
    nrcd=0

    for isf,batchNo in enumerate(batchList):
        # 140604
        if config.RunName.startswith('SAMPL4') or \
            config.RunName.startswith('focusedLib'):
            summPath = faahDir + 'ADV_summ_0000001.csv' 
        else:
            summPath = faahDir+('summ/%s_summ_%07d.csv' % (sys,batchNo))
            
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
                  
                if config.RunName.startswith('SAMPL4') or \
                    config.RunName.startswith('focusedLib'):
                    (ligRecept,ligand,e,eff,nvdv,ninter) = flds
                else:
                    # ['Exp32', '13651', 'xEyeSiteXtl5NI', 'ZINC00039702', '-4.63', '-0.421', '28', '4']
                    (expt,batch,ligRecept,ligand,e,eff,nvdv,ninter) = flds

                if ligRecept != exptRecept:
                    continue
                
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
            
            # Pass1: need to make this thresh criterion-sensitive, too
            
            if (dcrit == 'energy' and abs(dval) > TooBigE) or \
                (dcrit == 'ligEff' and abs(dval) > TooBigLigEff) or \
                (dcrit == 'ligEnth' and abs(dval) > TooBigLigEnth):
                noddE += 1
            else:
                
                # NB: dval used as generic value from here forward
                # "E" name is no longer correct, but doesn't hurt anything!
                
                allE.append(dval)
        inStr.close() # eo-batch file

    if len(allE)==0:
        print 'getThresh: No good summ files1?! Expt=%s NBatch=%d' % (exptName,len(batchList))
        return (0.,0.)
    
    # ASSUME all algorithms require n log n??
    allE.sort()
    if frac4Thresh==1.0:
        thresh = allE[-1]
    else:
        threshIdx = int(round(float(len(allE) * frac4Thresh)))
        thresh = allE[threshIdx]       
    
    if ncand > len(allE):
        print 'getThresh: ncand < allLig=%d; using all' % (len(allE))
        ncand = len(allE)-1
    bestCandThresh = allE[ncand]

    exptName = bldExptStr(exptKey)
    print 'getThresh: Expt=%s NRcd=%d NMissf=%d NOddE=%d ThreshFrac=%5.2f Thresh=%f BestCandThresh=%f' % \
            (exptName,nrcd,nmissf,noddE,frac4Thresh,thresh,bestCandThresh)
    
    return thresh, bestCandThresh
    
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

def ligIdx2zinc(ligIdx,zeroSuffix=False):
    '''NB: zeroIdx param OVERRIDES config.LigandZeroSuffix
    '''
    zid = ligIdx / 100
    sufId = ligIdx % 100

    if ligIdx < 0:
        if ligIdx in config.NZIdx2LigTbl:
            return config.NZIdx2LigTbl[ligIdx]
        else:
            print 'ligIdx2zinc: missing nonZincIdx?!',ligIdx
            return 'Lig??'

    # HACK for DUDE zincids (:
    # cf 1_AC_x1E66dude_ES_DD ZINC713_0 71300
    if config.RunName.startswith('DUDE') and zid < 800:
        zincid = 'ZINC%d' % (zid)
    else:
        zincid = 'ZINC%08d' % (zid)

    if sufId != 0 or not zeroSuffix:
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
    # DUDE accomodation
    
    if config.RunName.startswith('DUDE'):
        if lig.endswith('.VS'):
            lig = lig[:-3]
        elif lig.endswith('prasD'):
            lig = lig[:-5]
        elif lig.endswith('prasA'):
            lig = lig[:-5]

    if lig.find('ZINCCHEMBL') != -1:
        zincIdx = bldNonZincIdx(lig)
    elif lig.find('ZINC') != -1:
        
        # HACK: reduce ZINC id to integer, for reduced storage
        try:
            zincIdx = bldZincIdx(lig)
        except:
            print 'normLigand: odd ligand?!',lig
            # import pdb; pdb.set_trace()
            return -1
        
    else:
        zincIdx = bldNonZincIdx(lig)
        
    return zincIdx
    
def getInterDetails(itype,iinfo):
    ''' parse details of interactions from JSON dictionary
        cf. crawl_ADV.rptDataADV(), crawl_ADV.reducePlusInterDict()
        NB: PPI/TPI rcenter, ligcenter not included!
    '''
    if itype=='vdw':
        (rchain,raa,ratom) = iinfo
        liname = ''
    elif itype=='ppi' or itype=='tpi':
        (rchain,raa,rcenter,ligcenter) = iinfo
        liname = ''
        ratom = ''
    else:
        (rchain,raa,ratom,liname) = iinfo
        
    return (rchain,raa,ratom,liname)

def bldFocalInterTbl(faahDir,batch2ligTbl,activeLig=None):
    '''retrieval of interaction data for just focal ligands (collected by batch)
        NB: Exp47_LigandName_Hack removed
    '''
    
    # NB: receptInterTbl retains lig atom distinctions
    receptInterTbl = defaultdict( lambda: defaultdict(lambda: defaultdict( list ))) 
    # (rchain,raa,ratom) -> {itype - > {liname -> [ ligID ] } }
         
    nrcd = 0
    nmissf = 0
    nactive = 0
    
    for ib,bno in enumerate(batch2ligTbl.keys()):
        
        if config.RunName.startswith('SAMPL4') or \
            config.RunName.startswith('focusedLib'):
            inf = faahDir + 'ADV_inter_0000001.json' 
                           
        else:
            inf = faahDir+('inter/%s_inter_%07d.json' % (RunType,int(bno)))
            
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
            if config.RunName.startswith('SAMPL4') or \
                config.RunName.startswith('focusedLib'):
                (recept,ligand,interList) = interInfo
                expt = config.RunName
                ligIdx = normLigand(ligand)
            else:
                (expt,batchNo,recept,ligand, interList) = interInfo
            
                ligIdx = normLigand(ligand)                    
                if ligIdx not in focalLigTbl:
                    continue

            if activeLig != None:
                if ligIdx in activeLig:
                    nactive += 1
                            
            # don't check for SAMPL, as they use bldInterTbl_csv()

            for interTypeList in interList:
                itypeIdx, interEnum = interTypeList
                itype = InterTypes[itypeIdx]
                for iinfo in interEnum:
                    (rchain,raa,ratom,liname) = getInterDetails(itype,iinfo)
                            
                    k = (rchain,raa,ratom)
                    lk = liname
        
                    receptInterTbl[k][itype][lk].append( ligIdx )
                    
#                     if ligIdx<0 and itype != 'vdw':
#                         print '1,%s,%s,%s,%s,%s' % (k,itype,lk,ligand,ligIdx)

    if activeLig != None:
        print 'bldFocalInterTbl: NActive=%d' % (nactive)
                 
    print 'bldFocalInterTbl: NRcd=%d N(Chain+RAA+RAtom)=%d NMissFile=%d' % (nrcd,len(receptInterTbl),nmissf)

    ## NB: need to make serializable, for pickel!
    
    interTbl2 = {}
    for rfeat,ligTbl in receptInterTbl.items():         
        # (rchain,raa,ratom) -> {itype - > {latom -> [ ligID ] } }
        newLigTbl = {}
        for itype,ligFeatTbl in ligTbl.items():
            newLFTbl = {}
            for latom,ligIDList in ligFeatTbl.items():
                newLFTbl[latom] = ligIDList[:]
#                 for ligIdx in ligIDList:
#                     if ligIdx<0 and itype != 'vdw':
#                         print '2,%s,%s,%s,%s' % (rfeat,itype,latom,ligIdx)
                    
                
            newLigTbl[itype] = newLFTbl
            assert len(newLFTbl)==len(ligFeatTbl), 'bldFocalInterTbl: bad newLFTbl?!'
        assert len(newLigTbl)==len(ligTbl), 'bldFocalInterTbl: bad newLigTbl?!'
        interTbl2[rfeat] = newLigTbl
    assert len(interTbl2)==len(receptInterTbl), 'bldFocalInterTbl: bad interTbl2?!'
    
    return interTbl2

def analBldHIFeatures(faahDir,exptName,bnoList,thresh,featFile):    
    '''build energy discrimination predicate
        v3: use thresh, dont build explicit energy distribution
        use features built by '''
    
    print "analBldHIFeatures: %s NBatch=%d" % (exptName,len(bnoList))

    def tryExp47Hack(ligBits,oldLig,cc):
        try:
            newLig = ligBits[1]
        except Exception, e:
            print "Exp47?! %s oldLig='%s' ligBits='%s'" % (cc,oldLig,ligBits)
            newLig = None
        return newLig

    classTbl = {} # ligand -> aboveThreshP
    nrcd = 0
    noddE = 0
    npos = 0
    nmissf = 0
    
    # 2do HACK:  back-fitting SAMPLE directory structure
    if faahDir.find('sampl') != -1:
        SAMPL_expt = True
    else:
        SAMPL_expt = False
    
    ## first pass: build energy distribution
    
    for ib,bno in enumerate(bnoList):
        
        if SAMPL_expt:
            inf = faahDir+('summ/%s_summ_%05d.csv' % (RunType,bno))
        else:
            inf = faahDir+('summ/%s_summ_%07d.csv' % (RunType,bno))

        try:
            inStr = open(inf)
        except:
            nmissf += 1
            # import pdb; pdb.set_trace()
            continue
        
#         if (ib % 50) == 0:
#             print exptName,ib,inf,nrcd,len(classTbl)
        
        for il,line in enumerate(inStr.readlines()):
            if il == 0:
                continue
            nrcd += 1
            flds =line[:-1].split(',')

            try:
                if SAMPL_expt:
                    # ['s3VQ4_FBP', 's3VQ4', 'AVX101118_0', '-8.9', '-0.27', '57', '1']
                    (expt,recept,ligand,e,eff,nvdv,ninter) = flds
                else:
                    # ['Exp32', '13651', 'xEyeSiteXtl5NI', 'ZINC00039702', '-4.63', '-0.421', '28', '4']
                    (expt,batch,recept,ligand,e,eff,nvdv,ninter) = flds
            except Exception,e:
                print "analBldHIFeatures: bad line %s %d %s?!" % (inf,il,e)
                noddE += 1
                break # don't try to read other lines from inStr
                
                
            # Exp47!
            if Exp47_LigandName_Hack:
                ligBits = ligand.split('_')
                if ligBits[0].startswith('ZINC'):
                    ligand = ligBits[0]
                else:
                    callContext = ('analBldHIFeatures','%s,%d,%s' % (exptName,bno,il))
                    ligand = tryExp47Hack(ligBits,ligand,callContext)                    
                    if not ligand:
                        # errors counted in ??()
                        continue

            # HACK: reduce ZINC id to integer, for reduced storage
            
            if SAMPL_expt:
                zincIdx = ligand
            else:
                zincIdx = normLigand(ligand)
                        
            if abs(float(e)) > TooBigE:
                # print 'analFAAH: odd energy?!',isf,summF,il,e
                noddE += 1
                continue
                        
            if float(e)<=thresh:
                classTbl[zincIdx] = True
                npos += 1
                       
        inStr.close()

    print "analBldHIFeatures: %s NLig=%d NOddE=%d Thresh = %f NPos=%d NMissFile=%d" % (exptName,nrcd,noddE,thresh,npos,nmissf)    
        
    if nrcd==0:
        print 'analBldHIFeatures: no ligands?!'
        return None

    # (rchain,raa,ratom) -> {itype - > {latom -> [ ligID ] } }
    itpklf = InterTblDir + '%s_interTbl.pkl' % (exptName)
    interTbl = cPickle.load(open(itpklf,'rb'))                            

    featInfoTbl = {}
    ncnt = 0
    ndup = 0
    
    # 2do ASAP 141002: provide HIFLevel == 'LigSet' option;
    # merge all ligand atoms associated with same (rchain,raa,ratom,itype)
    # into rchain_raa_ratom_itype_AtomSetFeature
    
    # 140919
    # NB: multiple ratom, itype, latom for same rchain+raa can produce
    # DOUBLE COUNTS wrt/ same feature/ligand 
    # use separate fligTbl: f -> scoredLig table to ensure only first is counted
    fligTbl = {}

    for rfeat,ligTbl in interTbl.items():            
        # (rchain,raa,ratom) -> {itype - > {latom -> [ ligID ] } }
        (rchain,raa,ratom) = rfeat
                 
        for itype,ligFeatTbl in ligTbl.items():
                            
            for latom,ligIDList in ligFeatTbl.items():
                f = bldFeature2(rchain,raa,ratom,itype,latom)
                if f == '':
                    continue
                if f not in featInfoTbl:
                    featInfoTbl[f] = [0,0]
                for ligID in ligIDList:
                    
                    if f in fligTbl and ligID in fligTbl[f]:
                        ndup += 1
                        continue
                    
                    ncnt += 1
                    if f in fligTbl:
                        fligTbl[f][ligID] = True
                    else:
                        fligTbl[f] = {ligID: True}
                        
                    if ligID in classTbl:
                        featInfoTbl[f][0] += 1
                    else:
                        featInfoTbl[f][1] += 1

    print "analBldHIFeatures: %s NLig=%d NPos=%d NCnt=%d NDup=%d Nfeatures=%d" % \
        (exptName,nrcd,npos,ncnt,ndup,len(featInfoTbl))
        
    allFeat = featInfoTbl.keys()
    allFeat.sort()
    
    nbadEntropy = 0
    
    # Read by analCluster()
    outs = open(featFile,'w')
    outs.write('F,Pr,Info,PosEnt,NegEnt\n')
    for f in allFeat:
        # TP: in low-mode,     feature present
        # FN: in low-mode,     no feature
        # FP: not in low-mode, feature present
        # TN: not in low-mode, no feature
        tp = featInfoTbl[f][0]
        fn = npos - tp
        fp = featInfoTbl[f][1]
        tn = nrcd - npos - fp
        
        all = tp+fp+tn+fn
#         assert (fn>0 and tn>0), 'analBldHIFeatures: bad counts1?!'
#         assert all == nrcd, 'analBldHIFeatures: bad counts2?!'
        if not(fn>=0 and tn>=0 and all == nrcd):
            nbadEntropy += 1
            outs.write('%s,,,,,?1,%d,%d,%d,%d\n' % (f,tp,fp,npos,nrcd))
            continue

        try:
            ig, prFeat, posEnt,negEnt = infoGain(tp,fp,tn,fn)
        except Exception,ex:
            # print 'analHIFeatures: badEntropy?!',exptName,ex,ifeat,npos,nrcd,tp,fp,tn,fn
            nbadEntropy += 1
            outs.write('%s,,,,,?2,%d,%d,%d,%d\n' % (f,tp,fp,npos,nrcd))
            continue
                    
        outs.write('%s,%f,%f,%f,%f\n' % \
                   (f,prFeat,ig,posEnt,negEnt))
    outs.close()
    print "analBldHIFeatures: %s NBadEntropy=%d" % (exptName,nbadEntropy)

def norm3AA(aa3s):
    # normalize RAA: convert 3char AA to one char; use %03d for pos
    aa3 = aa3s[:3]
    if aa3 in AADict:
        aa1 = AADict[aa3]
    else:
        aa1 = 'X'
    try:
        iapos = int(aa3s[3:])
        if iapos >= 1000:
            iapos = 999
    except:
        iapos = 999
    aas = '%s%03d' % (aa1,iapos)
    return aas

def bldExptStr(exptKey):
    s = '_'.join(exptKey)
    return s 

def analAll_(exptTbl,faahDirParent,outDir,analFn,paramTbl):
    '''generic function iterating any analFn over all experiment batches
        paramTbl includes all cross-experiment constant params,
        augments these with all params named by exptFile column headers
    '''
    
    # 2do HACK:  back-fitting SAMPLE directory structure
    if faahDirParent.find('sampl') != -1:
        SAMPL_expt = True
    else:
        SAMPL_expt = False

    if 'exptSubset' in paramTbl:
        exptList = paramTbl['exptSubset']
    else:
        exptList = None
        
    allExpt = exptTbl.keys()
    allExpt.sort()
    for exptKey in allExpt:
        exptNo,prot,recept,site,lib = exptKey

        if exptList and not (('Exp'+exptNo) in exptList):
            print 'analAll_: skipping',exptKey
            continue
        
        if SAMPL_expt:
            faahDir = faahDirParent
        else:
            faahDir = faahDirParent + 'Exp%s/' % (exptNo)

        exptStr = bldExptStr(exptKey)

        exptData = exptTbl[exptKey]
        exptData.update(paramTbl)

        if 'thresh' in exptData:
            thresh = exptData['thresh']
        elif paramTbl['needThresh']:
            print 'analAll_: no thresh?!',exptStr
            continue

        batchList = ranges2list( [(exptData['bstart'],exptData['bend'])] )
        
        print 'analAll_: analyzing features %s %d batches "%s"' % \
            (exptStr,len(batchList),exptTbl[exptKey])
        
        analFn(faahDir,outDir,exptStr,batchList,exptData)
         
def analBestRLIF(exptTbl,faahDirParent,outDir,nbest=None,exptList=None):
    'Accumulate RLIF frequency stats for best ligands'

    allExpt = exptTbl.keys()
    allExpt.sort()
    
    for exptKey in allExpt:
        exptNo,prot,recept,site,lib = exptKey
        
        if exptList and exptNo not in exptList:
            continue
        
        if config.RunName.startswith('SAMPL4') or \
            config.RunName.startswith('focusedLib'):
            faahDir = faahDirParent
        else:
            faahDir = faahDirParent + 'Exp%s/' % (exptNo)

        exptData = exptTbl[exptKey]
        exptName = bldExptStr(exptKey)

        config.NonZincLigTbl = {}
        config.NZIdx2LigTbl = {}
        config.NNonZincLig = 0
                                
        lowef = LowEDir+ exptName + '_lowE.csv'

#        nactive = 0

        if nbest==None:
            ligTbl, batch2ligTbl = loadBestLig(lowef)
        else:
            ligTbl, batch2ligTbl = loadBestLig(lowef,maxLig=nbest)         

        itpklf = InterTblDir + '%s_interTbl.pkl' % (exptName)
        if os.path.isfile(itpklf):
            print 'analBestRLIF: interTbl exists, loading',itpklf
            interTbl = cPickle.load(open(itpklf,'rb'))
            
            nonzf = LowEDir + '%s_nonZincLig.csv' % (exptName)
            if os.path.exists(nonzf):
                reader = csv.DictReader(open(nonzf))
                for i,entry in enumerate(reader):
                    # Ligand,NZIdx
                    config.NonZincLigTbl[ entry['Ligand'] ] = int(entry['NZIdx'])
                    config.NZIdx2LigTbl [ int(entry['NZIdx']) ] =  entry['Ligand']
                config.NNonZincLig = len(config.NonZincLigTbl)
                print 'analBestRLIF: %d NonZincLig loaded from %s' % (config.NNonZincLig,nonzf)
            else:
                print 'analBestRLIF: no NonZincLig',exptName

        else:
            print 'analBestRLIF: interTbl does not exist, building',itpklf
            
#            interTbl = bldFocalInterTbl(faahDir,batch2ligTbl,activeLig)
            interTbl = bldFocalInterTbl(faahDir,batch2ligTbl)
            cPickle.dump(interTbl, open(itpklf,'wb'))
            
#             if NNonZincLig>0:
#                 # cf. bldNonZincIdx()
#                 nzligFile = LowEDir + '%s_nonZincLig.csv' % (exptName)
#                 print 'analBestRLIF: Writing %d NonZinc ligands to ' % (NNonZincLig,nzligFile)
#                 allLig = NonZincLigTbl.keys()
#                 allLig.sort()
#                 outs = open(nzligFile,'w')
#                 outs.write('Ligand,NZIdx\n')
#                 for lig in allLig:
#                     outs.write('%s,%d\n' % (lig,NonZincLigTbl[lig]))
#                 outs.close()
                
        # receptInterTbl[k][itype][lk].append( zincIdx )
        # (rchain,raa,ratom) -> {itype - > {liname -> [ ligID ] } }
        
        rlifTbl = defaultdict(int) # RLIF -> freq

        for rfeat,ligTbl in interTbl.items():            
            # (rchain,raa,ratom) -> {itype - > {latom -> [ ligID ] } }
            (rchain,raa,ratom) = rfeat
            for itype,ligFeatTbl in ligTbl.items():
                for latom,ligIDList in ligFeatTbl.items():
                    f = bldFeature2(rchain,raa,ratom,itype,latom)
                    if f.strip() == '':
                        continue
                    rlifTbl[f] += len(ligIDList)
        
        print 'analBestRLIF: %s NNonZinc=%d NRLIF=%d' % (exptName,config.NNonZincLig,len(rlifTbl))       
        
        allRLIF = rlifTbl.keys()
        allRLIF.sort(key= lambda f: rlifTbl[f],reverse=True)
        
        outf = RLIFDir + ('%s_rlif.csv' % exptName)
        outs = open(outf,'w')
        outs.write('RLIF,Freq\n')
        for f in allRLIF:
            outs.write('%s,%d\n' % (f,rlifTbl[f]))
        outs.close()

def loadBestLig(lowef,maxLig=-1,normalizeLigand=True):
    '''read energy, batch from lowef
    return ligTbl: ligIdx -> [e,batch] and batch2ligTbl: batch -> [ligs in batch]
    if maxLig != -1, only maxLig lowest energy ligands included
    150629: normLigand() applied
    '''

    ligTbl = {} # lig -> [e,batch]                                                                                              
    batch2ligTbl = defaultdict(list) # batch -> [lig]                                                                       

    reader = csv.DictReader(open(lowef))
    for i,entry in enumerate(reader):
        # cf analFAAH                                                                                                       
        #  fldNameLine = 'Batch,Ligand,E,Eff,Nvdw,Ninter\n'  
                                                                       
        # HACK: let batch tag along in ligTbl in first pass
        ligand = entry['Ligand']
        batchNo = entry['Batch']
        
        if normalizeLigand:
            ligIdx = normLigand(ligand)                                                     
            ligTbl[ ligIdx ] = [ float( entry['E'] ), batchNo ]
        else:
            ligTbl[ ligand ] = [ float( entry['E'] ), batchNo ]

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

    return ligTbl, batch2ligTbl
  
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
    
def bldFeature2(chain,raa,ratom,itype,ligatom):
    '''create features based on binary/wvdw, posOnly/not GLOBAL variables
       NB: may return "", eg, if binary and itype is vdw
       NB: unicode (from interact JSON) converted to str()
    '''

    f = ''
    
    if ADFeatures=='binary' and itype == 'vdw':
        return ''

    prefix = bldFeaturePrefix(str(chain),str(raa))

    if HIFLevel == 'PosOnlyHIF':          
        f = prefix
    elif HIFLevel == 'TypedPosHIF':
        f = prefix + '_' + str(itype)
    elif HIFLevel == 'ROnlyHIF':
        f = prefix + ('_'.join(['',str(ratom),str(itype)]) )
    else: # HIFLevel == 'Full':
        # 150929: drop ligand atom index
        if len(ligatom)==0:
            latype = ''
        else:
            latype = str(ligatom)[0]
        f = prefix + ('_'.join(['',str(ratom),str(itype),latype]) )

    return f

def feature2bits(f):
    'return [chain,Rpos,RAA,Ratom,IType,Latom];  AApos is int'
    bits = f.split('_')
    aa = bits[1]
    aapos = int(aa[:-1]) # drop 1-character residue
    aar = aa[-1]
    bits[1] = aapos
    bits.insert(2,aar)
    return bits
   
def loadSAMPLTrueTbl(inf):
    print 'loadSAMPLTrueTbl: loading',inf
    trueLigTbl = {} # ligIdx -> ligand name
    inStr = open(inf)
    for il,line in enumerate(inStr.readlines()):
        # AVX40944     AVX40944_0   C[N@@H+](CC1=CC=CC=C1)CC2=C(C(=C(C=C2)OC)CCC3=CC=CC(=C3)C(=O)[O-])C(=O)[O-]
        flds = line.split(' ')
        flds = [f for f in flds if f != ''] 
        ligand = flds[1]
        ligIdx = normLigand(ligand)
        trueLigTbl[ligIdx] = ligand
    inStr.close()
    print 'loadSAMPLTrueTbl: done. NTrueLig=%d' % len(trueLigTbl)
    return trueLigTbl # ligIdx -> ligand name

def analSAMPLTrueEnergy(faahDir,outDir,exptName,bnoList,paramTbl):    
    '''analyze SAMPL low energies with benefit of TRUE docking data
       analAll_() complient'''

    runType = paramTbl['runType']
    edecimal = paramTbl['edecimal']
    trueLigTbl = paramTbl['trueLigTbl']
    
    posEVec = []
    negEVec = []
    posEDistTbl = {}
    negEDistTbl = {}
    for ib,bno in enumerate(bnoList):
        inf = faahDir+'summ/%s_summ_%05d.csv' % (runType,bno)
        inStr = open(inf)
        
        for il,line in enumerate(inStr.readlines()):
            if il == 0:
                continue
            flds =line[:-1].split(',')
            # Expt,Recept,Ligand,E,Eff,Nvdw,Ninter
            (expt,recept,ligand,e,eff,nvdw,ninter) = flds
               
            e = float(e)         
            if abs(e) > TooBigE:
                # print 'analFAAH: odd energy?!',isf,summF,il,e
                continue

            rEstr = rndEstr(e,edecimal)
            
            if ligand in trueLigTbl:
                posEVec.append(e)
                if rEstr in posEDistTbl:
                    posEDistTbl[rEstr] += 1
                else:
                    posEDistTbl[rEstr] = 1         
            else:
                negEVec.append(e)
                if rEstr in negEDistTbl:
                    negEDistTbl[rEstr] += 1
                else:
                    negEDistTbl[rEstr] = 1
                    
        inStr.close()

    posAvg,posSD = basicStats(posEVec)
    negAvg,negSD = basicStats(negEVec)
    print "analSAMPLTrue: %s PosAvg=%f PosSD=%f NegAvg=%f NegSD=%f" % (exptName,posAvg,posSD,negAvg,negSD)
            
    plot.plot2ExptE(exptName,posEDistTbl,negEDistTbl,edecimal)

def analTrueEnergy(faahDir,outDir,exptName,ligTbl,activeIdxSet):    
    '''analyze low energies with benefit of knowledge of activeIdxSet
    NB: dropped edecimal rounding?'''
    
    posEVec = []
    negEVec = []
    posEDistTbl = defaultdict(int)
    negEDistTbl = defaultdict(int)
    
    for ligIdx in ligTbl.keys():
                       
            e = ligTbl[ligIdx][0]         
            if abs(e) > TooBigE:
                # print 'analFAAH: odd energy?!',isf,summF,il,e
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

def ligHIF2arff(faahDir,summDir,exptName,bnoList,paramTbl):
    '''create ARFF encoding of expt's ligands wrt/ HIF
       analAll_() complient
    '''

    runType = paramTbl['runType']
    trueLigTbl = paramTbl['trueLigTbl']
    featureSetName = paramTbl['featureSetName']
    infoThresh = paramTbl['infoThresh']
    # hifEncoder = paramTbl['hifEncoder']

    # 2do HACK:  back-fitting SAMPLE directory structure
    if faahDir.find('sampl') != -1:
        print 'ligHIF2arff: SAMPL_expt not supported(:'
        return None
        # SAMPL_expt = True
    else:
        SAMPL_expt = False

    outs = open(summDir+('arff/%s_%s.arff' % (exptName,featureSetName)),'w')
    
    ## get this experiment's HIF
    featureFile = HIFDir+('%s.csv' % exptName)
    hifTbl = {}
    reader = csv.DictReader(open(featureFile))
    # F,Pr,Info,PosEnt,NegEnt
    for i,entry in enumerate(reader):
        info = float(entry['Info'])
        if info > infoThresh:              
            hif = entry['F']
            hifTbl[hif] = info  
    
    print 'ligHIF2arff: NHIF=%d' % (len(hifTbl))
    # get ligands' energies
    ligETbl = {}
    for ib,bno in enumerate(bnoList):
        if SAMPL_expt:
            inf = faahDir+('summ/%s_summ_%05d.csv' % (runType,bno))
        else:
            inf = faahDir+('summ/%s_summ_%07d.csv' % (runType,bno))

        inStr = open(inf)
        
        for il,line in enumerate(inStr.readlines()):
            if il == 0:
                continue
            flds =line[:-1].split(',')
            if SAMPL_expt:
                (expt,recept,ligand,e,eff,nvdw,ninter) = flds
            else:
                # Expt,Batch,Recept,Ligand,E,Eff,Nvdw,Ninter
                (expt,batch,recept,ligand,e,eff,nvdw,ninter) = flds
               
            e = float(e)         
            if abs(e) <= TooBigE:
                ligETbl[ligand] = e             
        inStr.close()

    outs.write('@relation %s_%s\n' % (featureSetName,exptName))
    outs.write('@attribute e numeric\n')
    allHIF = hifTbl.keys()
    allHIF.sort()
    for hif in allHIF:
        outs.write('@attribute %s {True, False}\n' % hif)
    outs.write('@attribute class {True, False}\n')
    outs.write('@data\n')
    
    nrcd = 0
    for ib,bno in enumerate(bnoList):
        if not SAMPL_expt:
            inf = faahDir+('inter/%s_inter_%05d.json' % (runType,bno))
            
#         else:
#             inf = faahDir+('inter/%s_inter_%07d.json' % (runType,bno))

        try:
            inStr = open(inf)
            allInter = json.load(inStr)
            inStr.close()
        except Exception, e:
            print 'ligHIF2arff: bad inter JSON?!',exptName,bno,e
            continue

        currLigand = ''
        ligHIFTbl = {}

        for il,interInfo in allInter:
            nrcd += 1
            flds =line[:-1].split(',')
            
            # cf crawl_ADV.rptData_ADV()
            # [ [Expt,BatchNo,Recept,Lig, [IType,[InterEnum] ] ] ]
            (expt,recept,ligand, interList) = interInfo
            for interTypeList in interList:
                itypeIdx, interEnum = interTypeList
                itype = InterTypes[itypeIdx]
                for iinfo in interEnum:
                    (rchain,raa,ratom,liname) = getInterDetails(itype,iinfo)

            assert False, '2do: update ligLib2arff for new bldFeature2()'
            # need to pass adFeatures,posOnly into LigLib2arff
            # was: hif = bldFeature(rchain,raa,ratom,itype)
            # now: hif = bldFeature2(chain,raa,ratom,itype,ligatom,adFeatures,posOnlyHIF)
                    
            
            assert expt!='Exp47', 'bldInterTbl: retest Exp47_LigandName_Hack!'
#             if Exp47_LigandName_Hack:
#                 ligBits = ligand.split('_')
#                 if ligBits[0].startswith('ZINC'):
#                     ligand = ligBits[0]
#                 else:
#                     ligand = ligBits[1]

            if il==1:
                currLigand = ligand

            ligHIFTbl[hif] = True
        
        # ligStr = '%s,' % ligand
        ligStr = '%f,' % ligETbl[ligand]
        for hif in allHIF:
            ligStr += '%s,' % (hif in ligHIFTbl)
        pos = currLigand in trueLigTbl
        ligStr += '%s\n' % (pos)
        outs.write(ligStr)
            
    outs.close()


def loadFragClusters_v2(fragClustf): 
    '''load fragClustTbl: rlif -> frag -> centroid-frag
    HACK: 2 pass: first identify fragCenters, then map all others in clique to it
    '''

    frag2CtrTbl = {} # rlif -> frag -> fragCtr

    # pass1
    allCliques = {} # rlif -> cliqCtr
    cliqCtr = {} # cliqID -> fragCtr
    allFragC = []
    prevRLIF = None
    reader = csv.DictReader(open(fragClustf))
    for i,entry in enumerate(reader):
        # RLIF,CliqIdx,IsCtr,Frag
        rlif = entry['RLIF']
        cidx = entry['CliqIdx']
        ctr = entry['IsCenter']
        frag = entry['Frag']
        if rlif==prevRLIF:
            if ctr=='1':
                allFragC.append(frag)
                cliqCtr[cidx] = frag
        else:
            allCliques[prevRLIF] = (cliqCtr,allFragC)
            cliqCtr = {} # cliqID -> fragCtr
            allFragC = []
            prevRLIF = rlif
    allCliques[prevRLIF] = (cliqCtr,allFragC)   
         
    print 'loadFragClusters: NFragC=%d' % (len(allFragC))
    
    # pass2            
    frag2CtrTbl = defaultdict(dict) # rlif -> frag -> fragCtr
    reader = csv.DictReader(open(fragClustf))
    for i,entry in enumerate(reader):
        # Clique,IsCenter,Fragment
        rlif = entry['RLIF']
        cidx = entry['CliqIdx']
        frag = entry['Frag']
        frag2CtrTbl[rlif][frag] = cliqCtr[cidx]
            
    return frag2CtrTbl
    
def ligRFC2arff(exptName,activeIdxSet,ethresh=None):
    '''create ARFF encoding of expt's ligands wrt/ RLIF+CentroidFrag qualification
    ASSUME ligCoord already built, loaded from ligCoord.pkl
    use fragClustf to collapse fragments to cluster center
    rlifPartition==True when fragments have been clustered wrt/ RLIF
    output ARFF includes ligand names; filtered by weka
    '''

    R2FCPickleFile = R2FCDir + exptName + '_r2fc.pkl'
    print 'Loading r2fcPickle...',
    origFragClustTbl = cPickle.load(open(R2FCPickleFile,'rb'))
    # rlif -> cliqueIdx -> (ctrFrag, [cliqueFrags] )
    print 'done.'
    
    # Build rlif -> frag -> fragCtr
    fragClustTbl = {}
    for rlif in origFragClustTbl:
        fragClustTbl[rlif] = {}
        for cidx in origFragClustTbl[rlif]:
            ctr = origFragClustTbl[rlif][cidx][0]
            for frag in origFragClustTbl[rlif][cidx][1]:
                fragClustTbl[rlif][frag] = ctr
    
    lowef = LowEDir+ exptName + '_lowE.csv'
    # NB: no nbest passed to loadBestLig; ALL returned
    ligTbl, foo = loadBestLig(lowef)
    
    lcFile = LigCoordDir + '%s_ligCoord.pkl' % (exptName)
    print 'Loading ligCoordPickle...',lcFile
    ligCoordTbl = cPickle.load(open(lcFile,'rb'))
    
    print 'ligRFC2arff: NLig = %d' % (len(ligCoordTbl))
    if ethresh != None:
        print 'ligRFC2arff: ethresh=%f used' % (ethresh)
    allLig = ligCoordTbl.keys()
    allLig.sort()
    
    ## identify RAtom, RFQ features
    
    allRAtomsSet = set()  # allows enumeration when coding data

    # both ratomLigTbl and rfqTbl are dict for speedy ligand lookups
    ratomLigTbl = defaultdict(dict) 
    rfqTbl = defaultdict(dict) # ratom_frag_itype_latom -> {lig: True}
    nmissFrag = 0
    nhiELig = 0
    nlig = 0
    missFragTbl = defaultdict(int)
    for il,ligIdx in enumerate(allLig):
        
        if ethresh != None:
            e = ligTbl[ligIdx][0]
            if e > ethresh:
                nhiELig += 1
                continue
        nlig += 1
        allMatch = ligCoordTbl[ligIdx]
        allMatch.sort(key = lambda m: m['rlif'] )
        
        for mi,matchDict in enumerate(allMatch):
            rlif = matchDict['rlif']
            
            # bits = rlif.split('_')
            # ra = '_'.join(bits[:3])            
            # allRAtomsSet.add(ra)
            # ratomLigTbl[ra][ligIdx] = True

            # 150704: usef full rlif
            allRAtomsSet.add(rlif)
            ratomLigTbl[rlif][ligIdx] = True
                        
            frag = matchDict['frag']
            
            if rlif not in fragClustTbl or frag not in fragClustTbl[rlif]:
                # print 'ligRFC2arff: missing from fragClustTbl?!',il,lig,mi,rlif,frag
                nmissFrag += 1
                missFragTbl[rlif+'_'+frag] += 1
                continue
            else:
                cfrag = fragClustTbl[rlif][frag]                
            
            rfq = rlif + ('+%s' % cfrag)
            rfqTbl[rfq][ligIdx] = True
    
    if ethresh != None:
        print 'ligRFC2arff: ethresh=%f ==> nhiELig=%d NLig=%d' % (ethresh,nhiELig,nlig)
    print 'ligRFC2arff: NRA=%d NRFQ = %d NMissFrag=%d NUniqMiss=%d' % \
        (len(ratomLigTbl),len(rfqTbl),nmissFrag,len(missFragTbl))
    
    missFragFile = ArffDir+('%s_missFrag.csv' % (exptName))
    outs = open(missFragFile,'w')
    allMiss = missFragTbl.keys()
    allMiss.sort(key=lambda k: missFragTbl[k],reverse=True)
    outs.write('Miss,NMiss\n')
    for miss in allMiss:
        outs.write('%s,%d\n' % (miss,missFragTbl[miss]))
    outs.close()

    arrffFile = ArffDir+('%s.arff' % (exptName))
    outs = open(arrffFile,'w')
        
    outs.write('@relation %s\n' % (exptName))
    outs.write('@attribute ligand string\n')
    
    outs.write('@attribute e numeric\n')

    ## RAtom features
    
    prevRAtom = None
    allRAtomsList = list(allRAtomsSet)
    allRAtomsList.sort()
    for ra in allRAtomsList:
        outs.write('@attribute %s {0,1}\n' % ra)
       
    ## RFQ features      
    allRFQ = rfqTbl.keys()
    allRFQ.sort()

    nrfqFeature = 0
    nlofreqRFQ = 0
    allRFQList = []  # separate list, because some of rfq in allRFQ will be disqualified
    for rfq in allRFQ:
        bits = rfq.split('+')
        frag = bits[1]

        # Always drop low frequency fragments
        if len(rfqTbl[rfq]) < FragMinLigFreq:
            nlofreqRFQ += 1
            continue
        
        nrfqFeature += 1
        allRFQList.append(rfq)       
        outs.write('@attribute %s {0,1}\n' % rfq)

    print 'ligRFC2arff: NLoFreqRFQ=%d' % (nlofreqRFQ)
    
    outs.write('@attribute class {0,1}\n')

    outs.write('@data\n')
    
    nout = 0
    nactive = 0
    nhighE = 0
    nNoFeature = 0
#     for il,ligIdx in enumerate(allLigIdx):
#         zeroSuffix = (exptName.find('_PR_') == -1)
#         ligand = ligIdx2zinc(ligIdx,zeroSuffix)

    nbitsVec = []
    nhiELig = 0
    for il,ligIdx in enumerate(ligTbl):
        
#         if il % 1000 == 0:
#             print 'ligRFC2arff: writing data...',il

        e = ligTbl[ligIdx][0]
        if ethresh != None:
            if e > ethresh:
                nhiELig += 1
                continue
                    
        if ligIdx not in ligCoordTbl:
            nNoFeature += 1
            continue

                
        nout += 1
        ligStr = ''

        ligStr += '"%s",' % (ligIdx)
            
        ligStr += '%f,' % (e)
        
        nbitsSet = 0
        for ra in allRAtomsList:
            ligStr += '%d,' % (ligIdx in ratomLigTbl[ra])
            if ligIdx in ratomLigTbl[ra]:
                nbitsSet += 1
         
        for rfq in allRFQList:
            ligStr += '%d,' % (ligIdx in rfqTbl[rfq])
            if ligIdx in rfqTbl[rfq]:
                nbitsSet += 1
                
        nbitsVec.append(nbitsSet)
        
        # ligIdx = normLigand(ligand)
        # active = ligIdx in activeIdxSet
        active = ligIdx in activeIdxSet
        if active:
            nactive += 1

        ligStr += '%d\n' % (active)
        outs.write(ligStr)
            
    outs.close()        

    print 'ligRFC2arff: %s NHiELig=%d/%d NRFQ=%d NRAtom=%d NLow=%d NNoFeature=%d Nout=%d NActive=%d/%d' % \
        (exptName,nhiELig,len(ligTbl),len(rfqTbl),len(allRAtomsList),len(ligTbl),nNoFeature,nout,nactive,len(activeIdxSet))

    print 'ligRFC2arff: %s NRFQFeature=%d' % (exptName,len(allRFQList))
        
    avg,sd = basicStats(nbitsVec)
    print 'ligRFC2arff: NBitsSet Avg=%f SD=%f' % (avg,sd)   

def ligRFC2SpArff(exptName,origFragClustTbl,activeIdxSet,ethresh=None):
    '''create ARFF encoding of expt's ligands wrt/ RLIF+CentroidFrag qualification
    ASSUME ligCoord already built, loaded from ligCoord.pkl
    use fragClustf to collapse fragments to cluster center
    rlifPartition==True when fragments have been clustered wrt/ RLIF
    output ARFF includes ligand names; filtered by weka
    '''

#     R2FCPickleFile = R2FCDir + exptName + '_r2fc.pkl'
#     print 'Loading r2fcPickle...',
#     origFragClustTbl = cPickle.load(open(R2FCPickleFile,'rb'))
    # rlif -> cliqueIdx -> (ctrFrag, [cliqueFrags] )
    print 'done.'
    
    # Build rlif -> frag -> fragCtr
    fragClustTbl = {}
    for rlif in origFragClustTbl:
        fragClustTbl[rlif] = {}
        for cidx in origFragClustTbl[rlif]:
            ctr = origFragClustTbl[rlif][cidx][0]
            for frag in origFragClustTbl[rlif][cidx][1]:
                fragClustTbl[rlif][frag] = ctr
    
    lowef = LowEDir+ exptName + '_lowE.csv'
    # NB: no nbest passed to loadBestLig; ALL returned
    ligTbl, foo = loadBestLig(lowef)
    
    lcFile = LigCoordDir + '%s_ligCoord.pkl' % (exptName)
    print 'Loading ligCoordPickle...',lcFile
    ligCoordTbl = cPickle.load(open(lcFile,'rb'))
    
    print 'ligRFC2SpArff: NLig = %d' % (len(ligCoordTbl))
    if ethresh != None:
        print 'ligRFC2SpArff: ethresh=%f used' % (ethresh)
    allLig = ligCoordTbl.keys()
    allLig.sort()
    
    ## identify RAtom, RFQ features
    
    allRAtomsSet = set()  # allows enumeration when coding data

    # both ratomLigTbl and rfqTbl are dict for speedy ligand lookups
    ratomLigTbl = defaultdict(dict) 
    rfqTbl = defaultdict(dict) # ratom_frag_itype_latom -> {lig: True}
    nmissFrag = 0
    nhiELig = 0
    nlig = 0
    missFragTbl = defaultdict(int)
    for il,ligIdx in enumerate(allLig):
        
        if ethresh != None:
            e = ligTbl[ligIdx][0]
            if e > ethresh:
                nhiELig += 1
                continue
        nlig += 1
        allMatch = ligCoordTbl[ligIdx]
        allMatch.sort(key = lambda m: m['rlif'] )
        
        for mi,matchDict in enumerate(allMatch):
            rlif = matchDict['rlif']
            
            # bits = rlif.split('_')
            # ra = '_'.join(bits[:3])            
            # allRAtomsSet.add(ra)
            # ratomLigTbl[ra][ligIdx] = True

            # 150704: usef full rlif
            allRAtomsSet.add(rlif)
            ratomLigTbl[rlif][ligIdx] = True
                        
            frag = matchDict['frag']
            
            if rlif not in fragClustTbl or frag not in fragClustTbl[rlif]:
                # print 'ligRFC2arff: missing from fragClustTbl?!',il,lig,mi,rlif,frag
                nmissFrag += 1
                missFragTbl[rlif+'_'+frag] += 1
                continue
            else:
                cfrag = fragClustTbl[rlif][frag]                
            
            rfq = rlif + ('+%s' % cfrag)
            rfqTbl[rfq][ligIdx] = True
    
    if ethresh != None:
        print 'ligRFC2SpArff: ethresh=%f ==> nhiELig=%d NLig=%d' % (ethresh,nhiELig,nlig)
    print 'ligRFC2SpArff: NRA=%d NRFQ = %d NMissFrag=%d NUniqMiss=%d' % \
        (len(ratomLigTbl),len(rfqTbl),nmissFrag,len(missFragTbl))
    
    missFragFile = ArffDir+('%s_missFrag.csv' % (exptName))
    outs = open(missFragFile,'w')
    allMiss = missFragTbl.keys()
    allMiss.sort(key=lambda k: missFragTbl[k],reverse=True)
    outs.write('Miss,NMiss\n')
    for miss in allMiss:
        outs.write('%s,%d\n' % (miss,missFragTbl[miss]))
    outs.close()

    atIdx = {}
    arrffFile = ArffDir+('%s_sp.arff' % (exptName))
    outs = open(arrffFile,'w')
        
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

    outs.write('@attribute dummyVar string\n')
    atIdx['dummyVar'] = 0

    outs.write('@attribute ligand string\n')
    atIdx['ligand'] = 1
    
    outs.write('@attribute e numeric\n')
    atIdx['e'] = 2

    ## RAtom features
    
    prevRAtom = None
    allRAtomsList = list(allRAtomsSet)
    allRAtomsList.sort()
    
    begIdx = len(atIdx)
    for ra in allRAtomsList:
        outs.write('@attribute %s {0,1}\n' % ra)
        atIdx[ra] = begIdx
        begIdx += 1
        
    ## RFQ features      
    allRFQ = rfqTbl.keys()
    allRFQ.sort()

    nrfqFeature = 0
    nlofreqRFQ = 0
    allRFQList = []  # separate list, because some of rfq in allRFQ will be disqualified
    begIdx = len(atIdx)
    for rfq in allRFQ:
        bits = rfq.split('+')
        frag = bits[1]

        # Always drop low frequency fragments
        if len(rfqTbl[rfq]) < FragMinLigFreq:
            nlofreqRFQ += 1
            continue
        
        nrfqFeature += 1
        allRFQList.append(rfq)       
        outs.write('@attribute %s {0,1}\n' % rfq)
        atIdx[rfq] = begIdx
        begIdx += 1

    print 'ligRFC2SpArff: NLoFreqRFQ=%d' % (nlofreqRFQ)
    
    outs.write('@attribute class {0,1}\n')
    atIdx['class'] = len(atIdx)
    
    outs.write('@data\n')
    
    nout = 0
    nactive = 0
    nhighE = 0
    nNoFeature = 0
#     for il,ligIdx in enumerate(allLigIdx):
#         zeroSuffix = (exptName.find('_PR_') == -1)
#         ligand = ligIdx2zinc(ligIdx,zeroSuffix)

    nbitsVec = []
    nhiELig = 0
    for il,ligIdx in enumerate(ligTbl):
        
#         if il % 1000 == 0:
#             print 'ligRFC2arff: writing data...',il

        e = ligTbl[ligIdx][0]
        if ethresh != None:
            if e > ethresh:
                nhiELig += 1
                continue
                    
        if ligIdx not in ligCoordTbl:
            nNoFeature += 1
            continue

                
        nout += 1
        ligStr = '{'

        ligStr += '%d "%s",' % (atIdx['ligand'],ligIdx)
            
        ligStr += '%d %f,' % (atIdx['e'],e)
        
        nbitsSet = 0
        for ra in allRAtomsList:
            if ligIdx in ratomLigTbl[ra]:
                ligStr += '%d 1,' % (atIdx[ra])
                nbitsSet += 1
         
        for rfq in allRFQList:
            if ligIdx in rfqTbl[rfq]:
                ligStr += '%d 1,' % (atIdx[rfq])
                nbitsSet += 1
                
        nbitsVec.append(nbitsSet)
        
        # ligIdx = normLigand(ligand)
        # active = ligIdx in activeIdxSet
        active = ligIdx in activeIdxSet
        if active:
            nactive += 1

        ligStr += '%d %d}\n' % (atIdx['class'],active)
        outs.write(ligStr)
            
    outs.close()        

    print 'ligRFC2SpArff: %s NHiELig=%d/%d NRFQ=%d NRAtom=%d NLow=%d NNoFeature=%d Nout=%d NActive=%d/%d' % \
        (exptName,nhiELig,len(ligTbl),len(rfqTbl),len(allRAtomsList),len(ligTbl),nNoFeature,nout,nactive,len(activeIdxSet))

    print 'ligRFC2SpArff: %s NRFQFeature=%d' % (exptName,len(allRFQList))
        
    avg,sd = basicStats(nbitsVec)
    print 'ligRFC2SpArff: NBitsSet Avg=%f SD=%f' % (avg,sd)   
    
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

def bldRLIFTbl(rdf):
    
    ## first pass: build dict of ALL interaction features
    rlifInfoTbl = defaultdict(lambda: defaultdict(lambda: defaultdict(list))) # pos -> ratom -> itype -> [atoms] 
    
    reader = csv.DictReader(open(rdf))
    for i,entry in enumerate(reader):
        # cf analBestRLIF
        # RLIF,Freq
        rlif = entry['RLIF']  
        fbits = feature2bits(rlif)
        chain,rpos,raa,ratom,itype,latom = fbits
        
        if itype=='vdw' and not VDWIncuded:
            continue

        if chain not in ['A','B']:
            continue
        
        if chain == 'B' and AChainOnly:
            continue
        
#         if len(latom)>0:
#             latom = latom[0]
            
        rlifInfoTbl[rpos][ratom][itype].append(latom)

    ## second pass:             
    rlifTbl = {}
    for rpos in rlifInfoTbl:
        rposStr = '%03d' % (rpos)
        allRAA = rlifInfoTbl[rpos].keys()  # NB: required to make defaultdict keys stable!
        for raa in allRAA:
            allIT = rlifInfoTbl[rpos][raa].keys()
            for itype in allIT:
                f = ''
                if RLIFLevel=='PAI':
                    f = '_'.join([rposStr,raa,itype])
                    rlifTbl[f] = True
                elif RLIFLevel=='PAISet' or RLIFLevel=='PAIatom':
                    latomSet = rlifInfoTbl[rpos][raa][itype]
                    if RLIFLevel=='PAISet':
                        latomSet.sort()
                        atomSetString = '%s' % (latomSet)
                        f = '_'.join([rposStr,raa,itype,atomSetString])
                        rlifTbl[f] = True
                    else: # 'PAIatom'
                        for a in latomSet:
                            f = '_'.join([rposStr,raa,itype,a])
                            rlifTbl[f] = True

    return rlifTbl
      
def scoreLigList(ligList):
    '''compute unique, nactive over ligList
    '''
    
    ligSet = set(ligList)
    actSet = [lig for lig in ligSet if lig.find('CHEMBL') != -1]
    return len(ligSet), len(actSet) 
                                              
def bldRAtomFrag(ligCoordTbl,outf):
    '''report per-RAtom X itype sets of fragments, to be clustered
    '''
    
    # ligand -> [ {rlif,frag,fragIdx, latomFull,dropped,flaNames,fmap,laIdx,laCoord,fragCoord} ]
 
    print 'bldRAtomFrag: NLig = %d' % (len(ligCoordTbl))
    allLig = ligCoordTbl.keys()
    allLig.sort()

    railafTbl = defaultdict(lambda: defaultdict(list)) # rlif -> frag ->  [ligands]
    for lig in allLig:
        
        allMatch = ligCoordTbl[lig]
        allMatch.sort(key = lambda m: m['rlif'] )
        
        for mi,matchDict in enumerate(allMatch):
            rlif = matchDict['rlif']
            frag = matchDict['frag']
            railafTbl[rlif][frag].append(lig)
            
    print 'bldRAtomFrag: NRAIF = %d' % (len(railafTbl))
    
    outs = open(outf,'w')
    outs.write('RAI,Frag,NLig,NUniq,NActive\n')
    allRAILA = railafTbl.keys()
    allRAILA.sort()
    nout = 0
    fragFreqVec = []
    # print 'RAI,NFrag,NLig,NAct'
    for rai in allRAILA:
        fragFreqVec.append(len(railafTbl[rai]))
        totLig = 0
        totAct = 0
        for frag in railafTbl[rai]:
            ligList = railafTbl[rai][frag]
            nuniqLig,nactLig = scoreLigList(ligList)
            outs.write('%s,%s,%d,%d,%d\n' % (rai,frag,len(ligList),nuniqLig,nactLig))
            nout += 1
            totLig += len(railafTbl[rai][frag])
            totAct += nactLig
        # print '%s,%d,%d,%d' %  (rai,len(railafTbl[rai]),totLig,totAct)
    outs.close()
    avg,sd = basicStats(fragFreqVec)
    print 'bldRAtomFrag: NRAI=%d NRAIFOut=%d NFrag/RAI Avg=%f SD=%f' % (len(allRAILA),nout,avg,sd)

def bldRAtomFragSim(fragSimf,ratomFragf,ratomFSimf):
    '''combine RAI-qualified fragments with inter-frag similarites
    to create partition of inter-frag similarites  wrt/ RAI-qualified subsets
    used as basis for rai-frag clustering by bldFragClusters()
    '''
    
    fragSimTbl = defaultdict(dict)
    reader = csv.DictReader(open(fragSimf))        
    for i,entry in enumerate(reader):
        # Frag1,Frag2,Sim
        frag1 = entry['Frag1']
        frag2 = entry['Frag2']
        sim =   float(entry['Sim'])
        fragSimTbl[frag1][frag2] = sim
        
    raiFragSet = defaultdict(set) # rai -> [fragments]
    reader = csv.DictReader(open(ratomFragf))        
    for i,entry in enumerate(reader):
        # RAI,Frag,NLig,NUniq,NActive
        rai = entry['RAI']
        frag = entry['Frag']
        raiFragSet[rai].add(frag)
        
    outs = open(ratomFSimf,'w')
    outs.write('RAI,Frag1,Frag2,Sim')

    allRAI = raiFragSet.keys()
    allRAI.sort()
    for rai in allRAI:
        allFrag = list(raiFragSet[rai])
        allFrag.sort()
        for frag1 in allFrag:
            for frag2 in allFrag:
                if frag2 <= frag1:
                    continue
                try:
                    sim = fragSimTbl[frag1][frag2]
                except Exception,e1:
                    try:
                        sim = fragSimTbl[frag2][frag1]
                    except Exception,e2:
                        print 'huh?!',e2
                outs.write('%s,%s,%s,%f\n' % (rai,frag1,frag2,sim))
    
    outs.close()
          
    
# import psycopg2
        
           
def getAtoms(mol):
    # http://openbabel.org/docs/dev/UseTheLibrary/PythonDoc.html#using-iterators
    
    # Note that OBMolTorsionIter returns atom IDs which are off by one. 
    # That is, you need to add one to each ID to get the correct ID. 
    
    # cf http://forums.openbabel.org/Read-pdb-files-in-C-td4657503.html
    # http://forums.openbabel.org/Unexpected-behavior-with-GetResidue-td4657736.html
    
    atomList = []

    tbl = ob.OBElementTable()

#     atcntTbl = defaultdict(int)

    # print 'aidx,aindex,aid,atype,anum,aname1,aname2'

    for atom in ob.OBMolAtomIter(mol):
        aid = atom.GetId()
        atype = atom.GetType()
        r = atom.GetResidue()
        if r==None:
            aname = tbl.GetSymbol(atom.GetAtomicNum())
        else:
            aname = r.GetAtomID(atom).strip()

#         aidx = atom.GetIdx()
#         aindex = atom.GetIndex()
#         anum = atom.GetAtomicNum()
#         aname2 = tbl.GetSymbol(atom.GetAtomicNum())
            
#         print aidx,aindex,aid,atype,anum,aname1,aname2
            
        atomList.append( (int(aid), aname, atype) )
    

    return atomList        

def mapLig2Features(ligIdx,fragments,ligMol):

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
            
        fatoms = getAtoms(fragMol)
    
        # print 'FRAGMENT',fatoms
        bindDict[fi]['fragment'] = frag
        bindDict[fi]['fatoms'] = fatoms
    
        # TJO 11/6/2015
     '''new code block using smart pattern matching and mapping
        intended to produce same data structures and mapping, but
        with fewer (no!) missing mappings
        '''
    pat = ob.OBSmartsPattern()
    if  pat.Init(frag) and pat.Match(ligMol):

        #print 'NIso=%d' % (pat.NumMatches())
        if pat.NumMatches() == 0:
            # print "No isomorphs?!",zincid,fi,frag
            bindDict[fi]['maps'] = None
        else:

            maps = list()
            for p in pat.GetUMapList():
                maps.append( [(a,b-1) for (a,b) in enumerate(p)] )

            # print fi,ii, maps

            bindDict[fi]['maps'] = maps

        else:
            zincid = ligIdx2zinc(ligIdx)
            print "Error making mapper?!",zincid,fi,frag
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
            
        Subject:     Re: recapping - PS
        Date:     Wed, 25 Mar 2015 08:50:27 -0700
        From:     TJ ODonnell <tjo@acm.org>        
        
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
               
#     frag2check.sort(key = (lambda f: len(bindDict[f]['maps'])))
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
#             print '\n* incomplete assignment %s FragIdx=%d' % (zincid,f)
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
        
#         laIdxList = [lidx for fidx,lidx in useMap]
#         for (lidx,lname,ltype) in bindDict['latoms']:
#             if lidx in laIdxList:
#                 flaNames.append(lname)
            
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

def bldLig2frag(ligTbl,exptName,verbose=False):
    '''return ligFragTbl: ligIdx -> {fragment,mapIdx,nmiss,useMap,flaNames,dropped}
    doing RECAP fragmentation directly from PDBQT 
    vs. separate canon,RECAP steps as in _v1
    '''

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
    
    pat = ob.OBSmartsPattern();
   
    ngood=0
    npoor=0
    nslide=0
    nerr = 0
        
    for iz,ligIdx in enumerate(ligTbl.keys()):
                
        e,batch = ligTbl[ligIdx]
        # zincid (vs ligIdx) used for error messages
        zincid = ligIdx2zinc(ligIdx)
                
        if config.RunName.startswith('focusedLib'):
            ligand = ligIdx2zinc(ligIdx)
            pdbqf = ProcDir + LigPathTbl[ligand]
        else:
            bno = int(batch)
            pdbqf = dockFile(exptName,bno,ligIdx)
    
            if pdbqf == None:
                errs.write('missing pdbqt file: %d,%s\n' % (bno,ligIdx))
                dls = open(config.DropLigFile, 'a')
                dls.write('%s,bldLig2Frag: missing pdbqt file: %d,%s\n' % (zincid,bno,ligIdx))
                dls.close()
                nerr += 1
                continue
        
        ligMol = ob.OBMol()
        obc.ReadFile(ligMol,pdbqf)
        
        # NB RECAP works directly from ligmol    
        # don't need to build canonical string
        # except to include it in bindDict
        
        ligSmilesC = obc.WriteString(ligMol)
        # this returns both the canonSmiles string, but also PDBQT file name?!
        lsbits = ligSmilesC.split()
        canon = lsbits[0]

        # make copy in order to alter mol
        amol = ob.OBMol(ligMol)
        
        currRecap = recap.Recap(amol,4)
        for si,bondName in enumerate(currRecap.bondNames):
            pat.Init(currRecap.smarts[bondName])
            currRecap.apply(pat, si)
        currRecap.decide_multiples()
        currRecap.split()
        
        ligRecap = obc.WriteString(amol,True)
        # this returns both the canonSmiles string, but also PDBQT file name?!
        lrbits = ligRecap.split()
        recapStr = lrbits[0]
        
        fragments = recapStr.split('.')
           
        ###################################
        
        bindDict = mapLig2Features(ligIdx,fragments,ligMol)

        ###################################
        
        # NB: augment bindDict with canon here vs. in mapLig2Features, 
        # since obc created here
        
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
                               
    print 'bldLig2frag: NLig=%d NGood=%d NErr=%d NPoor=%d NSlide=%d' % \
        (len(ligTbl),ngood,nerr,npoor,nslide)
        
    outs.close()
    errs.close()
    
    return l2fTbl
   
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

def pybelBits2binary(fpbits):
    bitlist=list('0'*1024)
    for item in fpbits:
        bitlist[item-1]='1'
    bitstring = ''.join(bitlist)
    return bitstring

def bldfragFP(fragFile,fragFPFile):   
    
    outs = open(fragFPFile,'w')
    outs.write('Fragment,FP,FP2,FPnum\n')
    reader = csv.DictReader(open(fragFile))
    for i,entry in enumerate(reader):
        # fragment,nref
        frag = entry['fragment']
        mol = pybel.readstring('can',frag)
        fp = mol.calcfp()
        bits1 = pybelBits2binary(fp.bits)
        bits2 = ''.join([format(num,'032b') for num in fp.fp])
        numlist = [num for num in fp.fp]
        outs.write('%s,%s,%s,"%s"\n' % (frag,bits1,bits2,numlist))
    outs.close()

def rptFragSim(fragFile,fragSimFile,minSim=1e-3):
    '''filter from fragments to pairwise similarities
    reporting to file is important because of lig2fragments openbabel dependency
    '''
    
    reader = csv.DictReader(open(fragFile))
    molFPTbl = {}
    for i,entry in enumerate(reader):
        # fragment,nref
        frag = entry['Fragment']
        mol = pybel.readstring('can',frag)
        fp = mol.calcfp()
        molFPTbl[frag] = fp

    outs = open(fragSimFile,'w')
    outs.write('Frag1,Frag2,Sim\n')
    nsmall = 0
    nout = 0
    allFrag = molFPTbl.keys()
    allFrag.sort()
    for fi,frag1 in enumerate(allFrag):
        for fj,frag2 in enumerate(allFrag):
            if fj >= fi:
                continue
            fp1 = molFPTbl[frag1]
            fp2 = molFPTbl[frag2]
            tsim = fp1 | fp2
            if tsim < minSim:
                nsmall += 1
            else:
                outs.write('%s,%s,%f\n' % (frag1,frag2,tsim))
                nout += 1

    outs.close()
    print 'rptFragSim: NFrag=%d Nsmall (< %f) =%d NOut=%d' % (len(allFrag),minSim,nsmall,nout)
   
def bldLigCoordFromPDBQT(exptName,r2fTbl,ligTbl):
    ''' ligCoord : ligIdx -> [matchDict]
    each matchDict corresp to rlif+frag+ligand+fragIdx instance
    it extends r2fTbl's fields {rlif,frag,fragIdx,latomFull,dropped,flaNames,fmap}
    with {laIdx,laCoord,fragCoord=fragInfoList}
    '''

    ligAtoms = defaultdict(list) # lig -> [ {rlif,frag,fragIdx,latomFull,dropped,flaNames,fmap} ]
    for rlif in r2fTbl.keys():  # RLIF -> frag -> [ (zincid,fragIdx,latomFull,fragInfo) ]
        for frag in r2fTbl[rlif]:
            for ligFragInfo in r2fTbl[rlif][frag]:
                zincid,fragIdx,latomFull,fragInfo = ligFragInfo
                
                ligAtoms[zincid].append( {  'rlif': rlif,
                                            'frag': frag,
                                            'fragIdx': fragIdx,
                                            'latomFull': latomFull,
                                            # NB: map,dropped all the way from assignFrag()
                                            #     then via lig2frag and rlif2frag!
                                            'dropped': fragInfo['dropped'],
                                            'flaNames': fragInfo['flaNames'],
                                            'fmap': fragInfo['useMap']} )
                
    print 'bldLigCoordFromPDBQT: %s NLigands=%d' % (exptName,len(ligAtoms))                           
    
    allLig = ligAtoms.keys()
    allLig.sort()
    
    errf = LigCoordDir + exptName + 'ligCoord_err.txt'
    errs = open(errf,'w')
    nerr = 0
    
    ligCoord = defaultdict(list)
    # maintain freq counts of number of RLIF-interacting latoms associated with same ligand
    multMatchFreq = defaultdict(int)
    nmissla = 0
    for il,zincid in enumerate(allLig): 

        # each match corresp to rlif+frag+ligand+fragIdx instance
        multMatchFreq[len(ligAtoms[zincid])] += 1
#             for matchDict in ligAtoms[lig]:
#                 print 'bldLigCoordFromPDBQT:',exptName,lig,matchDict['rlif'],matchDict['fragIdx'],matchDict['latomFull'],matchDict['frag']
            
        ligIdx = normLigand(zincid)
                    
        if ligIdx not in ligTbl:
            errs.write('ligand missing from lowLig %s %s\n' % (exptName,zincid))
            nerr += 1
            continue
        
        e,batch = ligTbl[ligIdx]
                
        if config.RunName.startswith('focusedLib'):
            pdbqf = ProcDir + LigPathTbl[zincid]
        else:
            bno = int(batch)
            pdbqf = dockFile(exptName,bno,ligIdx)
    
            if pdbqf == None:
                errs.write('ligand missing pdbqt file %s %s %s\n' % (exptName,bno,zincid))
                nerr += 1
                continue

        allMol = pybel.readfile('pdbqt', pdbqf)
        # print 'len(allMol)',len(allMol)
        ligMol = allMol.next() # ASSUME only one mol PDBQT
        obmol = ligMol.OBMol
        
        for matchDict in ligAtoms[zincid]:
            laIdx = None
            laCoord = None
            laName = matchDict['latomFull']
            # Leave room for any atoms dropped in assignFrag()
            nfragAtoms = len(matchDict['fmap'])
            fragInfoList = [None for i in range(nfragAtoms) ] # (atomName,idx,coord)
            l2fDict = {}
            for fidx,lidx in matchDict['fmap']:
                if lidx == None:  # dropped
                    continue
                l2fDict[lidx] = fidx
            
            # 2do: any way to replace iteration thru all atoms with index addressing?!
            for res in ob.OBResidueIter(obmol):
                for obatom in ob.OBResidueAtomIter(res):
                    pbatom = pybel.Atom(obatom)
                    idx = pbatom.idx
                    # print pbatom.idx,atName,pbatom.coords,pbatom.type
                    if idx-1 in l2fDict:
                        # NB: need to strip GetAtomID(obatom) !
                        atName = res.GetAtomID(obatom).strip()

                        fragInfoList[ l2fDict[idx-1] ] = (atName,idx,pbatom.coords)
                            
                        if atName == laName:
                            laIdx = idx
                            laCoord = pbatom.coords

            if laIdx == None:
                nmissla += 1
                errs.write('cant find ligAtom %s %s %s\n' % (exptName,zincid,matchDict))
                nerr += 1
                continue
            
            matchDict['laIdx'] = laIdx
            matchDict['laCoord'] = laCoord
            matchDict['fragCoord'] = fragInfoList
            
            for fi,finfo in enumerate(fragInfoList):
                if finfo == None:
                    fragInfoList[fi] =  (matchDict['flaNames'][fi],None,None)
                else:
                    # assert finfo[0] == matchDict['flaNames'][fi]
                    if finfo[0] != matchDict['flaNames'][fi]:
                        print 'huh'
            
            ligCoord[ligIdx].append(matchDict)
    
    errs.close()
    
    print 'bldLigCoordFromPDBQT: %s NMissLigAtom=%d NErr=%d' % (exptName,nmissla,nerr)
    
#     print 'bldLigCoordFromPDBQT: Number of RLIF-interacting latoms associated with same ligand'
#     allFreq = multMatchFreq.keys()
#     allFreq.sort()
#     print 'NLAtom,NLig'
#     for f in allFreq:
#         print '%d,%d' % (f,multMatchFreq[f])
    print 'bldLigCoordFromPDBQT: multMatchFreq: "%s"' % (str(multMatchFreq))

    ## need real dictionary to serialize/pickel!
    ligCoord2 = {} 
    for ligIdx in ligCoord.keys():
        ligCoord2[ligIdx] = ligCoord[ligIdx][:]
    return ligCoord2

def loadActives(actives_file):
    actives = []
    fs=open(actives_file,'r')
    for active in fs.readlines():
        actives.append(active.strip())
    fs.close()
    return actives
            
### from bind2RLIF2frag

def dockFile(exptName,bno,ligIdx):
    # .../Dock/1_RT_x2ZD1_RT_NNRTI_NNRTInADJ_NNRTI_DD/0000010_ZINC06556034_0_out_Vina_VS.pdbqt


    # DUDE hack!
    # cf normLigand()
    if config.RunName.startswith('DUDE'):
        if exptName.find('_PR_') != -1:  # all DUDE but PR use zeroSuffix
            zsuf = True
        else:
            zsuf = False
            
        ligand = ligIdx2zinc(ligIdx,zeroSuffix=zsuf)
    
        if exptName.find('_PR_') != -1:  # all DUDE but PR use zeroSuffix
            if ligand.find('CHEMBL') != -1:
                ligand += 'prasA'
            else:
                ligand += 'prasD'
    else:
        ligand = ligIdx2zinc(ligIdx)
            
    if config.RunName.startswith('IN_LEDGF_47-51'):
    # 150607:  2do: getPDBQT is adding cruff ?!
    #         .../Dock/47_IN_x3NF6_A_IN_LEDGF_LEDGF_NF/0046745_fahv.x3NF6_A_IN_LEDGF_ZINC17212212_1071250201_out_Vina_VS.pdbqt
    #         .../Dock/47_IN_x3NF6_A_IN_LEDGF_LEDGF_NF/0046745_                      ZINC17212212_           out_Vina_VS.pdbqt

        filePat = ProcDir + '/%s/%07d_fahv.*_%s_*_out_Vina_VS.pdbqt' % (exptName,bno,ligand)
        mfiles = glob.glob(filePat)

        if len(mfiles) == 0:
            print 'bldRLIF2frag: missing PDBQT?!',ligand
            return None
        elif len(mfiles)==1:
            pdbqf = mfiles[0]
   
    elif config.RunName.startswith('SAMPL4'):
        # 150613
        # .../processed/sampl/2_In_s3ZCM_A_LEDGF_SAMPL/0000714_AVX38779_2_out_Vina_VS.pdbqt
        # .../processed/sampl/Y3/s3NF8_B/AVX38779_2_1/AVX38779_2_1_VINoutput_Vina_VS.pdbqt

        # ebits = faahAnal.exptName2bits(exptName)
        ebits = exptName.split('_')
        # 1_In_s3NF8_A_LEDGF_SAMPL
        exptNo,protein,recept,chain,site,SAMPL = ebits
        recept = recept + '_' + chain
        dockDir2 = ProcDir + '%s/%s/%s/' % (site,recept,ligand)
        fname = ligand + '_VINoutput_Vina_VS.pdbqt'     
        pdbqf =  dockDir2 + fname

    elif config.RunName.startswith('DUDE'):
        dockDir2 =   ProcDir + '%s/' % (exptName)
        
        if exptName.find('_AC_') != -1 or \
           exptName.find('_AD_') != -1 or \
           exptName.find('_AM_') != -1 or \
           exptName.find('_HM_') != -1 or \
           exptName.find('_PG_') != -1:
            
            fname = ('%07d_' % bno) + ligand + '.VS.pdbqt'

        else:
            fname = ('%07d_' % bno) + ligand + '_out_Vina_VS.pdbqt'
            
        pdbqf =  dockDir2 + fname

    else:
        fname = ('%07d_' % bno) + ligand + '_out_Vina_VS.pdbqt'     
        pdbqf =  ProcDir + fname
    
    if os.path.exists(pdbqf):
        return pdbqf
    else:
        return None
    
def rptUniqLigFrag(l2fFile,ufragFile,minRef=1):

    reader = csv.DictReader(open(l2fFile))
    fragTbl = defaultdict(list) # frag -> [ligand]
    totFrag = 0
    for i,entry in enumerate(reader):
        # Ligand,FragIdx,Fragment,Map
        ligand = entry['Ligand']
        frag = entry['Fragment']
        fragTbl[frag].append(ligand)
        totFrag += 1
    
    allFrag = fragTbl.keys()
    
    nout = 0
    ndrop = 0
    allFrag.sort() # sort lexically on fragment
    outs = open(ufragFile,'w')
    outs.write('Fragment\n')
    for frag in allFrag:
        if len(fragTbl[frag]) > minRef:
            outs.write('%s\n' % (frag))
            nout += 1
        else:
            ndrop += 1

    outs.close()
    print 'rptUniqLigFrag: NFrag=%d TotFragRef=%d NLowFreq=%d NOut=%d' % \
        (len(fragTbl),totFrag,ndrop,nout)
        
    return nout
 
def rptFragLig(l2fFile,ufragFile,minRef=1):

    reader = csv.DictReader(open(l2fFile))
    fragTbl = defaultdict(list) # frag -> [ligand]
    totFrag = 0
    for i,entry in enumerate(reader):
        # Ligand,FragIdx,Fragment,Map
        ligand = entry['Ligand']
        frag = entry['Fragment']
        fragTbl[frag].append(ligand)
        totFrag += 1
    
    allFrag = fragTbl.keys()
    
    nout = 0
    allFrag.sort() # sort lexically on fragment
    outs = open(ufragFile,'w')
    outs.write('Fragment,Ligand,\n')
    for frag in allFrag:
        if len(fragTbl[frag]) > minRef:
            nout += 1
            for ligand in fragTbl[frag]:
                outs.write('%s,%s\n' % (frag,ligand))

    outs.close()
    print 'rptUniqLigFrag: NFrag=%d TotFragRef=%d NOut=%d' % \
        (len(fragTbl),totFrag,nout)
    
def rptR2FTbl(r2fTbl,outf):
    
    allRLIF = r2fTbl.keys()  # RLIF -> frag -> [ (zincid,fragIdx) ]
    allRLIF.sort()
    outs = open(outf,'w')
    outs.write('RLIF,Fragment,NLig,"[(ZINCID,FragIdx)]"\n')
    for rlif in allRLIF:
        for frag in r2fTbl[rlif]:
            outs.write('%s,%s,%d,"%s"\n' % (rlif,frag,len(r2fTbl[rlif][frag]),r2fTbl[rlif][frag]))
    outs.close()

def addAtomIdx(pdbqf):

    allMol = pybel.readfile('pdbqt', pdbqf)
    # print 'len(allMol)',len(allMol)
    ligMol = allMol.next() # ASSUME only one mol PDBQT
    obmol = ligMol.OBMol
    
    # 2do: any way to replace iteration thru all atoms with index addressing?!
    for res in ob.OBResidueIter(obmol):
        for obatom in ob.OBResidueAtomIter(res):
            pbatom = pybel.Atom(obatom)
            idx = pbatom.idx
            # print pbatom.idx,atName,pbatom.coords,pbatom.type
            # NB: need to strip GetAtomID(obatom) !
            atName = res.GetAtomID(obatom).strip()
            laCoord = pbatom.coords

                    
def bldRLIF2frag(exptName,ligTbl,l2fTbl,errf,dockFileTbl=None,bstart=None,bend=None,verbose=False):
    '''Aggregate RECAP fragments associated with RLIF
    
    RLIF -> frag -> [ (zincid,fragIdx,map,latomFull) ]
    (Re-) acquire RLIF directly from PDBQT using crawlADV routines
    
    dockFileTbl added to make combineLigLibTarget() work
    '''

    nmissFrag = 0
    nmissflig = 0
    nflig = 0
    r2fTbl = defaultdict( lambda: defaultdict(list) ) # RLIF -> frag -> [ (zincid,fragIdx,latomFull,fragInfo) ]

    ligList = ligTbl.keys()
    ligList.sort()
    
    nmissIdx = 0
    
    errs = open(errf,'w')
    nerr = 0
    
    for iz,ligIdx in enumerate(ligList):
        zincid = ligIdx2zinc(ligIdx)
        e,batch = ligTbl[ligIdx]
        
        if zincid not in l2fTbl:
            # print 'bldRLIF2frag: covTbl lig missing from fragments?!',lig
            nmissflig += 1
            continue
        
        nflig += 1
        ligFrags = l2fTbl[zincid]

        if dockFileTbl:
            pdbqf = ProcDir + dockFileTbl[zincid]
            
        elif config.RunName.startswith('focusedLib'):
            pdbqf = ProcDir + LigPathTbl[zincid]
            
        elif config.RunName.startswith('SAMPL4'):
            # NB: no need for batchno in SAMPL4
            foo=999
            pdbqf = dockFile(exptName,foo,ligIdx)
            if pdbqf == None:
                errs.write('SAMPL4 cant find PDBQT %s\n' % (zincid))
                nerr += 1
                continue

        else:
            bno = int(batch)
            pdbqf = dockFile(exptName,bno,ligIdx)
    
            if pdbqf == None:
                errs.write('missing pdbqt file %s %s %s\n' % (exptName,bno,ligIdx))
                nerr += 1
                continue
                    
        ligData = crawl_ADV.parseADPDBQT_ADV(pdbqf)
        
        # ligData:  {'recept': 'x3KFN_prASw0c0', 
        #            'hba': [('B', 'ILE50', 'N', 'O4'), ('B', 'ASP29', 'N', 'O6'), ('B', 'ASP30', 'N', 'O6'), ('A', 'ASP29', 'N', 'O2')],
        #            'hbd': [('A', 'GLY27', 'O', 'N1')],
        #            'vdw': [('B', 'VAL82', 'CG1') ...] }
        #            'e': -13.3, 
        #            'leff': -0.369    
        #            'nresult' : 1    
        #            'recept': 'x3KFN_prASw0c0'  
        #            'src' : 's> 9 '               

        # 2do: tpi, ppi binary interactions refer to ligand "CENTER" coordinates      
        for itype in ['hba','hbd','mtl']:
            
            if itype not in ligData:
                continue
            
            for ir,rawRlif in enumerate(ligData[itype]):
                
                chain,res,ratom,latomFull = rawRlif
                       
                # 150516: drop atomIdx!
                latom = latomFull[0]
                # 2do: replace with ASSERT
                try:
                    latomFoo = int(latomFull[1]) 
                except Exception,e:
                    nmissIdx += 1
                    errs.write('odd latom %s %s %s\n' % (zincid,ir,rawRlif))
                    nerr += 1
                    continue
                
                rlif = bldFeature2(chain,res,ratom,itype,latom)
                
                fragFnd = False
                for fi,fragInfo in enumerate(ligFrags):
                    if fragInfo == None:
                        errs.write('no fragInfo %s %s\n' % (zincid,fi))
                        nerr += 1
                        continue
                    if fragInfo['flaNames'] == None:
                        errs.write('no flaNames %s %s %s\n' % (zincid,fi,fragInfo))
                        nerr += 1
                        continue
                    
                    ## FOUND: (full ligand atom name from docking found in 
                    ##    fragInfo['flaNames'] associated with ligand's fragments by bldLig2Frag-> mapLig2Features -> getAtoms() ->
                    ##        OBRESIDUE.GetAtomID(atom).strip() OR OBElementTable.GetSymbol(atom.GetAtomicNum()) 
                    ## all fragInfo from l2fTbl carried forward
                    if latomFull in fragInfo['flaNames']:
                        frag = fragInfo['fragment']
                        r2fTbl[rlif][frag].append( (zincid,fi,latomFull,fragInfo) )
                        fragFnd = True
                        break
                    ## 
                    
                if not fragFnd:
                    errs.write('latomFull not found %s %s\n' % (zincid,latomFull))
                    nerr += 1
                    nmissFrag += 1
                    continue
                
        if verbose and iz % 1000 == 0:
            print 'bldRLIF2frag: iz=%d nunboundLigAtom=%d' % (iz,nmissFrag)
            
        ## eo-ligList
        
    errs.close()
    
    print 'bldRLIF2frag: NLig=%d NMissLig=%d NMissFrag=%d NErr=%d NMissIdx(SAMPL4)=%d' % (nflig,nmissflig,nmissFrag,nerr,nmissIdx)
    ## need real dictionary to serialize/pickel!
    
    r2fTbl2 = {} # RLIF -> frag -> [ (zincid,fragIdx,latomFull,fragInfo) ]
    for rlif in r2fTbl.keys():
        newFragTbl = {}
        for frag in r2fTbl[rlif].keys():
            newFragTbl[frag] = r2fTbl[rlif][frag][:]
        r2fTbl2[rlif] = newFragTbl

    return r2fTbl2

def bldR2FCtr_v2(exptName,r2fTbl,ligTbl,outf,verbose=False):
    ''' returns: r2fcTbl: rlif -> cliqueIdx -> (ctrFrag, [cliqueFrags] )
    with rlif-specific fragment clustering
    symmetric entries both kept for bldFragClusters()
    outputs fragClustf with all RLIFs clusters
    '''

    minNFragLig = 2
    
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
    outs.write('Expt,RLIF,NFrag,NSimFrag,NClust\n')
    # r2fTbl: RLIF -> frag -> [ (zincid,fragIdx,map,latomFull) ]
    allRLIF = r2fTbl.keys()
    allRLIF.sort()
    r2fcTbl = {}  # rlif -> cliqueIdx -> (ctrFrag, [cliqueFrags] )
    if verbose:
        print 'bldR2FCtr_v2-verbose: NRLIF=%d' % (len(allRLIF))
        print 'Idx,RLIF,NFrag,Sec'
    for ir,rlif in enumerate(allRLIF):
        allFrag = r2fTbl[rlif].keys()
        allFrag.sort()
        nlowFreq = 0

        if verbose:
            begTime = datetime.datetime.now()
        
        fragSimTbl =  defaultdict(dict) # frag1 -> frag2 -> tsim
        for frag1 in allFrag:

            nfragLig = len(r2fTbl[rlif][frag1])
            if nfragLig < minNFragLig:
                nlowFreq += 1
                continue
            
            mol1 = pybel.readstring('can',frag1)
            fp1 = mol1.calcfp()
            for frag2 in allFrag:
                if frag2 <= frag1:
                    continue
                
                nfragLig = len(r2fTbl[rlif][frag2])
                if nfragLig < minNFragLig:
                    continue

                mol2 = pybel.readstring('can',frag2)
                fp2 = mol2.calcfp()
                tsim = fp1 | fp2
                
                fragSimTbl[frag1][frag2] = tsim
                # NB: symmetric entries both kept for bldFragClusters()
                fragSimTbl[frag2][frag1] = tsim
             
        if len(fragSimTbl)==0:
            r2fcTbl[rlif] = {}
            outs.write('%s,%s,%d,%d,%d\n' % (exptName,rlif,len(allFrag),0,0))
            continue 
            
        
        rlifClustTbl = bldFragClusters(fragSimTbl,simThresh=simThresh,distThresh=distThresh)
        # rlifClustTbl: cliqueIdx -> (ctrFrag, [cliqueFrags] )
        
        outs.write('%s,%s,%d,%d,%d\n' % (exptName,rlif,len(allFrag),len(fragSimTbl),len(rlifClustTbl)))

        r2fcTbl[rlif] = rlifClustTbl
        
        if verbose:
            elapTime = datetime.datetime.now() - begTime
            print '%d,%s,%d,%s' % (ir,rlif,len(allFrag),elapTime.seconds)
        

        
    outs.close()
    return r2fcTbl

def loadFragDist(inf,verboseFreq=None):
    ''' build fragDist: frag1,frag2 -> (1.0 - Tanimoto similarity)
    
    loadFragDist built from output of query:
    select a.fragment afrag, b.fragment bfrag ,oc.tanimoto(a.fp,b.fp)
    from rik.fragdudefp2 a, rik.fragdudefp2 b
    where oc.nbits_set(a.fp) > 0 and oc.nbits_set(b.fp) > 0 and a.fragment
    
    very large for DUDE_PR; required lig2fragment.rptCullSimHist()
    
    '''
        
    fragDist = {}
    reader = csv.DictReader(open(inf))
    for i,entry in enumerate(reader):
        if verboseFreq != None and ((i % verboseFreq) == 0):
            print 'loadFragDist: Line:', i
            
        # generated by lig2fragment.rptCullSimHist()
        # NB: retaining \t separator used by pg_dump
        # 'Frag1\tFrag2\tSim\n' 
        frag1 = entry['Frag1']
        frag2 = entry['Frag2']
        
        # fragSim[frag1][frag2] = float(entry['Sim'])
        dist = 1 - float(entry['Sim'])
        if frag1 in fragDist:
            fragDist[frag1][frag2] = dist
        else:
            fragDist[frag1] = { frag2: dist}
                
        
    print 'loadFragDist: NFrag=%d' % (len(fragDist))
    return fragDist
 
def loadFragSim(inf):
    ''' build fragSim: frag1,frag2 -> Tanimoto similarity
    
    (n choose 2) upper diagonal pairs only, sorted on frag
    built via:
            
    select fragment, ob.fp(fragment) into fragdudefp2  from rik.fragdude ;
    
    select a.fragment, b.fragment, oc.tanimoto(a.fp,b.fp)
    from rik.fragdudefp2 a, rik.fragdudefp2 b 
    where (oc.nbits_set(a.fp & ~b.fp) + 
           oc.nbits_set(b.fp & ~a.fp) +
           oc.nbits_set(a.fp & b.fp)) > 0;
       
    meaning ALL pairs listed!
    '''
    
    ## first pass: get all fragments
    fragKeys = defaultdict(int)
    reader = csv.DictReader(open(inf))
    for i,entry in enumerate(reader):
        # fragment1,fragment2,tanimoto
        frag1 = entry['fragment1']
        frag2 = entry['fragment2']
        if frag1==frag2:
            continue
        fragKeys[frag1] += 1
        fragKeys[frag2] += 1
        
    allFrag = fragKeys.keys()
    allFrag.sort()
    
    # 2do: change to test/assert
    
    # NB: two comparisons not including self
    for frag in allFrag:
        if fragKeys[frag] != 2 * (len(allFrag) -1):
            print 'loadFragSim: odd number of comparisons?!',frag
    
    ## 2d pass: only include upper diagonal
    fragSim = defaultdict( lambda: defaultdict(float) )
    reader = csv.DictReader(open(inf))
    for i,entry in enumerate(reader):
        # fragment1,fragment2,tanimoto
        frag1 = entry['fragment1']
        frag2 = entry['fragment2']
     
        fidx1 = allFrag.index(frag1)
        fidx2 = allFrag.index(frag2)
        # only include upper diagonal
        if fidx1 > fidx2:
            continue
        
        fragSim[frag1][frag2] = float(entry['tanimoto'])
        
    print 'loadFragSim: NFrag=%d' % (len(allFrag))
    return fragSim
 
def bldfragDistArr(fragDistDict):
    '''produce numpy array of fragment distances
    '''

    allFrag = fragDistDict.keys()
    allFrag.sort()
    nfrag = len(allFrag)

    ## NB: allLig and allRLIF sorted above!
    spCov = dok_matrix((nfrag,nfrag), dtype=np.float)
    for fi,frag1 in enumerate(allFrag):
        for fj,frag2 in enumerate(fragDistDict[frag1].keys()):
                spCov[fi,fj] = fragDistDict[frag1][frag2]

    nnz1 = spCov.nnz
    if nnz1==0:
        print 'bldfragArray: nnz1==zero!?'
        return None

    # print 'bldfragArray: NNZ1=%d converting dok -> csr...' % (nnz1)
    spCov = spCov.tocsr()
    nnz2 = spCov.nnz
    if nnz2==0:
        print 'bldfragArray: nnz2==zero!?'
        return None

    # print 'bldfragArray: NNZ2=%d making dense...' % (nnz2)
    denseCov = spCov.toarray()
    nnz3 = np.count_nonzero(denseCov)
    if nnz3==0:
        print 'bldfragArray: nnz3==zero!?'
        return None
    # print 'bldfragArray: NNZ3=%d' % (nnz3)

    return denseCov
    
def splitRecap(allRecapf,canonDir):        
    '''utility to split allRecapf into individual experiment recap files 
    after processing via database recap
    '''
    
    assert False, '2do NEXT!'
    ins = open(allRecapf)
    prevExpt = None
    outs = None
    # 2do HACK: don't bother to break lines up into fields with csv.DictReader
    for il,line in enumerate(ins.readlines()):
        # Expt,ZINCID,Batch,E,CanonSmiles
        cpos = line.find(',')
        exptName = line[:cpos]
        if prevExpt:
            outs.close()
            outs = open(canonDir+exptName+'_recap.csv','w')
            outs.write('ZINCID,Batch,E,CanonSmiles\n')


def bldFragClusterFile(fragSimf,outf,simThresh=0.5,distThresh=0.5):
    '''cluster fragments using fastcluster and DISTANCE flattening
    built from FILE
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
                sim = fragTbl[frag1][frag2]
                tot += sim
            if tot > maxSim:
                maxSim = tot
                bestFrag = frag1
            tot2 += tot
            
        if bestFrag==None:
            print 'getCenter: all sim==0.?!'
            # arbitrarily pick first frag
            bestFrag = fragList[0]
                
#         print 'getCenter: NFrag=%d NSim=%d bestAvgSim=%f AvgSim=%f' % \
#             (nfrag,nsim,maxSim/(nfrag-1),tot2/nsim)

        return bestFrag

    fragTbl = defaultdict( lambda: defaultdict(float) ) # frag1 -> frag2 -> tsim
    fragSet = set()
    reader = csv.DictReader(open(fragSimf))
    
    for i,entry in enumerate(reader):
        # Frag1,Frag2,Sim
        frag1 = entry['Frag1']
        frag2 = entry['Frag2']
        fragSet.add(frag1)
        fragSet.add(frag2)
        if frag1==frag2:
            continue
        sim = float(entry['Sim'])
        if sim > simThresh:
            fragTbl[frag1][frag2]= sim
            # NB: symmetric entries in fragTbl used by getCenter()
            fragTbl[frag2][frag1]= sim
            
        

    # http://danifold.net/fastcluster.html
    # The argument X is either a compressed distance matrix or a collection of 
    # N observation vectors in D dimensions as an (N x D) array. 
    # Apart from the argument preserve_input, the methods have the same input 
    # and output as the functions of the same name in the package scipy.cluster.hierarchy. 
    # Therefore, I do not duplicate the documentation and refer to the SciPy documentation 
    # http://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html and 
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html for 
    # further details.   
          
    allFrag = list(fragSet)
    nfrag = len(allFrag)
    allFrag.sort()
    fragIdxTbl = {}

    fragDistArr = np.ones((nfrag,nfrag))
    for fi, frag in enumerate(allFrag):
        fragIdxTbl[frag] = fi
        fragDistArr[fi,fi] = 0.
    
    nnonz = 0
    
    for frag1 in fragTbl.keys():
        fidx1 = fragIdxTbl[frag1]
        # self-distance = zero
        for frag2 in fragTbl[frag1].keys():
            fidx2 = fragIdxTbl[frag2]
            sim = fragTbl[frag1][frag2]
            # NB: convert to DISTANCE
            dist = 1. - sim
            fragDistArr[fidx1,fidx2] = dist
            # NB: need to redundantly put both symmetric distances or squareform() bitches
            fragDistArr[fidx2,fidx1] = dist
            nnonz += 1
    sparsity = nnonz / (nfrag*(nfrag-1)/2.)
    nnonz2 = np.count_nonzero(fragDistArr)

    distVec = distance.squareform(fragDistArr)
    maxDist = max(distVec)
    
    print 'bldFragClusterFile: NFragment = %d NNonZero=%d (sparsity=%f) maxDist=%f' % (nfrag,nnonz,sparsity,maxDist)

    ## using sklearn KMeans
#     km = sklearn.cluster.KMeans(n_clusters=nclust)
#     km.fit(fragArr)
#     labels = km.labels_
#     centroids = km.cluster_centers_
    
    clustering = fastcluster.linkage(distVec,method='ward')
    # NB: sch.linkage doesn't support Ward method when "raw" observations not available
    # clustering = sch.linkage(distVec,method='average')
    
#     flatClust_n = sch.fcluster(clustering, nclust, 'maxclust')
#     flatClust_d = sch.fcluster(clustering, 0.5*maxDist, 'distance')
#     flatClust_i = sch.fcluster(clustering, 0.5, 'inconsistent')
#     flatClust = flatClust_n

#     if flatIdx == 'nclust':
#         flatClustOpt = [nclust, 'maxclust']
#     elif flatIdx == 'dist':
#         flatClustOpt = [0.1, 'distance'] # simThresh/4
#     elif flatIdx == 'incon':
#         flatClustOpt = [0.5, 'inconsistent']
#           
#     flatClust = sch.fcluster(clustering, flatClustOpt[0], flatClustOpt[1])

    flatClust = sch.fcluster(clustering, distThresh, 'distance')
    
    # clustRoot,clustTree = sch.to_tree(clustering, rd=True)
    
    # 2do:  better way to get most-central fragment for each cluster?
    # http://stackoverflow.com/questions/9362304/how-to-get-centroids-from-scipys-hierarchical-agglomerative-clustering
    
    clustTbl = defaultdict(list) # clustID -> [fragments]
    for fragIdx,cidx in enumerate(flatClust):
        clustTbl[cidx].append(allFrag[fragIdx])
    
    clustIndices = clustTbl.keys()
    print 'bldFragClusterFile: NClust found=%d' % (len(clustIndices))
    
    cliqueSizes = [len(clustTbl[c]) for c in clustTbl]
    print 'bldFragClusterFile: NClique=%d %s' % (len(clustTbl),cliqueSizes)

    clustIndices.sort()
    centralFrag = {} # clusterIdx -> most central frag
    for cidx in clustIndices:
        ctrFrag = getCenter(clustTbl[cidx])
        centralFrag[cidx] = ctrFrag    

    outs = open(outf,'w')
    outs.write('Clique,IsCenter,Fragment\n')
                                 
    for cidx in clustIndices:
        for frag in clustTbl[cidx]:
            isCenter = (frag==centralFrag[cidx])
            outs.write('%d,%d,%s\n' % (cidx,isCenter,frag))    
    outs.close() 

def bldFragClusters(fragSimTbl,simThresh=0.7,distThresh=0.5):
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
                
#         print 'getCenter: NFrag=%d NSim=%d bestAvgSim=%f AvgSim=%f' % \
#             (nfrag,nsim,maxSim/(nfrag-1),tot2/nsim)

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
    
    ## using sklearn KMeans
#     km = sklearn.cluster.KMeans(n_clusters=nclust)
#     km.fit(fragArr)
#     labels = km.labels_
#     centroids = km.cluster_centers_
    
    try:
        clustering = fastcluster.linkage(distVec,method='ward')
    except Exception,e:
        print e
        import pdb; pdb.set_trace()
        
    # NB: sch.linkage doesn't support Ward method when "raw" observations not available
    # clustering = sch.linkage(distVec,method='average')
    
#     flatClust_n = sch.fcluster(clustering, nclust, 'maxclust')
#     flatClust_d = sch.fcluster(clustering, 0.5*maxDist, 'distance')
#     flatClust_i = sch.fcluster(clustering, 0.5, 'inconsistent')
#     flatClust = flatClust_n

#     if flatIdx == 'nclust':
#         flatClustOpt = [nclust, 'maxclust']
#     elif flatIdx == 'dist':
#         flatClustOpt = [0.1, 'distance'] # simThresh/4
#     elif flatIdx == 'incon':
#         flatClustOpt = [0.5, 'inconsistent']
#           
#     flatClust = sch.fcluster(clustering, flatClustOpt[0], flatClustOpt[1])

    flatClust = sch.fcluster(clustering, distThresh, 'distance')
    
    # clustRoot,clustTree = sch.to_tree(clustering, rd=True)
    
    # 2do:  better way to get most-central fragment for each cluster?
    # http://stackoverflow.com/questions/9362304/how-to-get-centroids-from-scipys-hierarchical-agglomerative-clustering
    
    clustTbl = defaultdict(list) # clustID -> [fragments]
    for fragIdx,cidx in enumerate(flatClust):
        clustTbl[cidx].append(allFrag[fragIdx])
    
    clustIndices = clustTbl.keys()
    
#     cliqueSizes = [len(clustTbl[c]) for c in clustTbl]
#     print 'bldFragClusters: NClique=%d %s' % (len(clustTbl),cliqueSizes)

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

def tstDockFile(ligTbl,exptName,verbose=False):
    '''ala bldLig2frag(), but just make sure dockFile() works!
    '''

    nerr = 0
    for iz,ligIdx in enumerate(ligTbl.keys()):
                
        e,batch = ligTbl[ligIdx]
                
        if config.RunName.startswith('focusedLib'):
            ligand = ligIdx2zinc(ligIdx)
            pdbqf = ProcDir + LigPathTbl[ligand]
        else:
            bno = int(batch)
            pdbqf = dockFile(exptName,bno,ligIdx)
        if pdbqf == None:
            nerr += 1
            
    print 'tstDock: %s Nerr=%d/%d' % (exptName,nerr,len(ligTbl))

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

        import pybel
        ob = pybel.ob
        import recap3 as recap
        import fastcluster  # only when using bldFragClusters()

    elif HostName.startswith('mjq'):
        print 'running local on mjq'
        BaseDir = '/home/Data/coevol-HIV/WCG/'

        import pybel
        ob = pybel.ob

    else:
        print 
        sys.exit( ('unknown host %s' % (HostName)) )

    config.RunName = 'DUDE_151005' # 'IN-LEDGF' 'SAMPL4' 'DUDE' 'focusedLib' 'SAMPL4' 

    ProcDir = BaseDir + 'processed/%s/'  % (config.RunName)
    
    CrawlDir = BaseDir + 'crawl/'
        
    SummRptDir = BaseDir + 'anal/%s/'  % (config.RunName)

    RunType = 'ADV'
    FeatureType = 'RLIF' # 'sampl' 'HIF'
    ADFeatures = 'binary' # 'wvdw' 'binary'
    HIFLevel =  'Full' # vs 'PosOnlyHIF' vs. 'ROnlyHIF' vs 'Full'

    Exp47_LigandName_Hack = False

    ExptFile = SummRptDir+'exptBatch.csv'

    FocalExpts = None 

    VDWIncuded=False
    AChainOnly=False
    RLIFLevel = 'PAIatom'  # 'PAI' vs. 'PAIset' vs 'PAIatom'

    HIFLevel =  'Full' # 'PosOnlyHIF' vs. ROnlyHIF' vs 'Full'

    if FeatureType == 'RLIF':
        FeatLbl = 'R_%s_%s' % (ADFeatures[0],RLIFLevel[3])
    elif FeatureType == 'HIF':
        FeatLbl = 'H_%s_%s' % (ADFeatures[0],HIFLevel[0],HIFLevel[0])
    elif FeatureType == 'sampl':
        FeatLbl = '%s_%s_%s' % (FeatureType[0],ADFeatures[0].upper(),HIFLevel[0])
    else:
        sys.exit( ('unknown FeatureType?! %s' % (FeatureType)) )
        
    NBestLig = None
    LigPathTbl = None

    simThresh = 0.5

    nclust = 50
    
    PlotFmt = 'pdf'

    RestrictFragCount = False
    MaxFrag = 1000
    config.LigandZeroSuffix=False
      
    ############################################################
    ## config.RunName-specific variations
    
    if config.RunName.startswith('DUDE'):
        
        ProcDir = BaseDir + 'processed/DUDE/'
        if HostName.startswith('hancock') or HostName.startswith('mjq'):
            CrawlDir = BaseDir + 'crawl/DUDE/'
        
        ActiveDir = SummRptDir+('actives/')
        ActiveNameTbl = {'1_AC_x1E66dude_ES_DD': 'aces',
                        '1_AD_x2E1Wdude_AX_DD': 'ada',
                        '1_AM_x1L2Sdude_PC_DD': 'ampc',
                        '1_CD_x1H00dude_K2_DD': 'cdk2',
                        '1_ED_x2RGPdude_FR_DD': 'egfr',
                        '1_HM_x3CCWdude_DH_DD': 'hmdh',
                        '1_PG_x3LN1dude_H2_DD': 'pgh2',
                        '1_PR_x3KF0_prASw0c0_AS_DD': 'pras',
                        '1_PR_x3KF0_prASw1c0_AS_DD': 'pras',
                        '1_PY_x1D3Gdude_RD_DD': 'pyrd',
                        '1_RO_x2ETRdude_CK_DD': 'rock1',
                        '1_RT_x2ZD1_RT_NNRTI_NNRTInADJ_NNRTI_DD': 'rtnn'
                        }
        
        # HACK: DUDE defaults to LigandZeroSuffix = True, but PR_3KF0 experiments set it back!
        config.LigandZeroSuffix = True
        # NBestLig = 1000

    if config.RunName.startswith('SAMPL4'):
        FocalExpts = ['1','2','3','4','5','6']  # LEDGF site only
 
    ## focusedLib depends on dockings' path files written by crawl_adv.Local_focusedLib_Top_visit_ADV()
    if config.RunName.startswith('focusedLib'):
        ligPathFile = CrawlDir+'filePaths.csv'
        LigPathTbl = loadLigPaths(ligPathFile)
        
        # NBestLig = 1000
        
    ############################################################
    
    print '<analFAAH config.RunName=%s begTime=%s>' % (config.RunName,datetime.datetime.now().strftime('%y%m%d_%H%M%S'))
    print '\thost=%s\n\tWCGData=%s\n\tsummRpt=%s\n' % (HostName, CrawlDir, SummRptDir)
    print '\texptFile=%s\n\texptList=%s\n\trunType=%s\n\tfeatureType=%s' % (ExptFile,FocalExpts,RunType,FeatureType)
    print '\tADFeatures=%s\n\tHIFLevel=%s\n\tFeatLbl=%s\n' % (ADFeatures,HIFLevel,FeatLbl)
    print '\tNBestLig=%s\n' % (NBestLig)

    ############################################################

    ## Top-level directories; created as needed below as part of each phase
    
    LowEDir = SummRptDir+'lowE/'
    InterTblDir = SummRptDir+'InterTable/'
    RLIFDir = SummRptDir+'RLIF_F/'
    L2FDir = SummRptDir+ 'L2F/'
    FragSimDir = SummRptDir+ 'FragSim/'
    R2FDir = SummRptDir + 'R2F/'
    LigCoordDir = SummRptDir+'LigCoord/'
    R2FCDir = SummRptDir + 'R2FC/'
    PlotDir = SummRptDir + 'plots/'
    ArffDir = SummRptDir+'ARFF/'
    HIFDir = SummRptDir+'HIF/'


# #     ############################################################
     
    print 'analFAAH: *PHASE 0: building exptTbl'
          
    exptTblPkl0 = SummRptDir+'exptTblU_v0.pkl'
    if os.path.exists(exptTblPkl0):
        print '## Loading pre-pickled exptTbl_v0',exptTblPkl0
        exptTbl0 = cPickle.load(open(exptTblPkl0,'rb'))
    else:
        print '## <bldExptTbl>'
        exptTbl0 = bldExptTbl(ExptFile)
        print '## </bldExptTbl>'
        print '## Saving pickled newExptTbl',exptTblPkl0
        cPickle.dump(exptTbl0, open(exptTblPkl0,'wb'))
      
    if FocalExpts != None:
        print '## Focused experiments only:',FocalExpts
      
    if not os.path.isdir(LowEDir):
        print 'analFAAH: creating lowE directory',LowEDir
        os.makedirs(LowEDir)
      
    exptTblPkl1 = SummRptDir+'exptTblU_v1.pkl'
    if os.path.exists(exptTblPkl1):
        print '## Loading pre-pickled exptTbl_v1',exptTblPkl1
        exptTbl = cPickle.load(open(exptTblPkl1,'rb'))
    else:
        print '## <analFAAH_expt>'
        exptTbl = analFAAH_expt(exptTbl0,FocalExpts,CrawlDir,SummRptDir,frac4Thresh=1.0,ncand=1000)
        print '## </analFAAH_expt>'
        print '## Saving pickled exptTbl_v1',exptTblPkl1
        cPickle.dump(exptTbl, open(exptTblPkl1,'wb'))
      
    # sys.exit('done building exptTbl')
    
#     ############################################################
     
    print 'analFAAH: *PHASE 1: building RLIF'
           
    exptTblPkl1 = SummRptDir+'exptTblU_v1.pkl'
    exptTbl = cPickle.load(open(exptTblPkl1,'rb'))
         
    if not os.path.isdir(InterTblDir):
        print 'analFAAH_1: creating InterTbl directory',InterTblDir
        os.makedirs(InterTblDir)
       
    if not os.path.isdir(RLIFDir):
        # NB: presence of RLIFDir taken as indication this phase accomplished!
        print 'analFAAH_1: creating RLIF directory',RLIFDir
        os.makedirs(RLIFDir)
       
        print '## <analBestRLIF>'
        if NBestLig == None:
            analBestRLIF(exptTbl,CrawlDir,SummRptDir,exptList=FocalExpts)
        else:
            analBestRLIF(exptTbl,CrawlDir,SummRptDir,exptList=FocalExpts,nbest=NBestLig)
        print '## </analBestRLIF>'
            
    else:
        print 'analFAAH_1: RLIF already constructed'
       
#     sys.exit('done building RLIF')
    
#     ############################################################
    # NB: subsequent phases all assume allExpt
        
    allExpt = exptTbl.keys()
    allExpt.sort()
              
    ############################################################
      
    print 'analFAAH: *PHASE 1-1: compute HIF'
    if not os.path.isdir(HIFDir):
        print 'analFAAH_1-1: creating HIFDir directory',HIFDir
        os.makedirs(HIFDir)

    for expt in allExpt:
        
        exptNo,prot,recept,site,lib = expt
        exptName = bldExptStr(expt)
        
        hifFile = HIFDir+('%s.csv' % (exptName))

        if os.path.exists(hifFile):
            print 'analFAAH_1-1: %s HIF file exist; skipping' % (exptName)
            continue

        if config.RunName.startswith('SAMPL4') or \
            config.RunName.startswith('focusedLib'):
            faahDir = CrawlDir
        else:
            faahDir = CrawlDir + 'Exp%s/' % (exptNo)

        batchList = ranges2list( [(exptTbl[expt]['bstart'],exptTbl[expt]['bend'])] )
        frac4Thresh = 0.02
        ncand = 1000
        dcrit='energy'
        
        thresh, bestCandThresh = getThresh(expt,RunType,batchList,faahDir,frac4Thresh,ncand,dcrit)
        
        analBldHIFeatures(faahDir,exptName,batchList,thresh,hifFile)

    # sys.exit('done building HIF')
        
    ############################################################
      
    print 'analFAAH: *PHASE 2: compute lig2frag'
    if not os.path.isdir(L2FDir):
        print 'analFAAH_2: creating L2FDir directory',L2FDir
        os.makedirs(L2FDir)
       
    for expt in allExpt:
        exptNo,prot,recept,site,lib = expt
        exptName = bldExptStr(expt)
  
        L2FPickleFile =   L2FDir + exptName + '_lig2Frag.pkl'
        if os.path.exists(L2FPickleFile):
            print 'analFAAH_2: %s L2F pkl exist; skipping' % (exptName)
            continue
                       
        begTime = datetime.datetime.now()
        begTimeStr = begTime.strftime('%y%m%d_%H%M%S')
  
        print '<analFAAH_2 %s %s>' % (exptName,begTimeStr)
            
        ## NonZincLigTbl is experiment-specific
        ## needs to be kept for reloading along with pickeled
        config.NonZincLigTbl = {} # ligand name -> nonzincID 
        config.NZIdx2LigTbl = {} # nonzincID -> ligand name
        config.NNonZincLig = 0
        if config.RunName.startswith('DUDE'):
            nonzf = LowEDir + '%s_nonZincLig.csv' % (exptName)
        else:
            nonzf = LowEDir + 'nonZincLig.csv'
        if os.path.exists(nonzf):
            # cf. analFAAH_expt()
            reader = csv.DictReader(open(nonzf))
            for i,entry in enumerate(reader):
                # Ligand,NZIdx
                config.NonZincLigTbl[ entry['Ligand'] ] = int(entry['NZIdx'])
                config.NZIdx2LigTbl [ int(entry['NZIdx']) ] =  entry['Ligand']
            config.NNonZincLig = len(config.NonZincLigTbl)
            print 'analFAAH_2: %d NonZincLig loaded from %s' % (config.NNonZincLig,nonzf)
        else:
            print 'analFAAH_2: no NonZincLig',exptName
    
        lowef = LowEDir+ exptName + '_lowE.csv'
        # NB: full set of ligands need to be fragmented for classification testing
        if NBestLig == None or (config.RunName.startswith('DUDE') or config.RunName.startswith('SAMPL4')):
            ligTbl, batch2ligTbl = loadBestLig(lowef)
        else:
            ligTbl, batch2ligTbl = loadBestLig(lowef,maxLig=NBestLig)
                                
        BindPPFile = L2FDir + exptName + '_lig2fragPP.txt'
    
        # tstDock(ligTbl,exptName)
        
        l2fTbl = bldLig2frag(ligTbl,exptName)
               
        print 'L2FPickle dumping to',L2FPickleFile
        cPickle.dump(l2fTbl, open(L2FPickleFile,'wb'))
                
    # sys.exit( ('Done compute lig2frag' ) )
                  
    ############################################################
  
    print 'analFAAH: *Phase 3: build bldRLIF2frag'
    
    if not os.path.isdir(R2FDir):
        print 'analFAAH_3: creating RLIF2F directory',R2FDir
        os.makedirs(R2FDir)
           
    # nbest = 1000
    for expt in allExpt:  
        exptNo,prot,recept,site,lib = expt
        exptName = bldExptStr(expt)
  
        R2FPickleFile =   R2FDir + exptName + '_rlif2frag.pkl'
  
        if os.path.exists(R2FPickleFile):
            print 'R2FPickleFile exists; using',R2FPickleFile
            continue    
       
        begTime = datetime.datetime.now()
        begTimeStr = begTime.strftime('%y%m%d_%H%M%S')
        print '<analFAAH_3  %s %s>' % (exptName,begTimeStr)
                
        ## NonZincLigTbl is experiment-specific
        ## needs to be kept for reloading along with pickeled
        config.NonZincLigTbl = {} # ligand name -> nonzincID 
        config.NZIdx2LigTbl = {} # nonzincID -> ligand name
        config.NNonZincLig = 0
        if config.RunName.startswith('DUDE'):
            nonzf = LowEDir + '%s_nonZincLig.csv' % (exptName)
        else:
            nonzf = LowEDir + 'nonZincLig.csv'
        if os.path.exists(nonzf):
            # cf. analFAAH_expt()
            reader = csv.DictReader(open(nonzf))
            for i,entry in enumerate(reader):
                # Ligand,NZIdx
                config.NonZincLigTbl[ entry['Ligand'] ] = int(entry['NZIdx'])
                config.NZIdx2LigTbl [ int(entry['NZIdx']) ] =  entry['Ligand']
            config.NNonZincLig = len(config.NonZincLigTbl)
            print 'analFAAH_3: %d NonZincLig loaded from %s' % (config.NNonZincLig,nonzf)
        else:
            print 'analFAAH_3: no NonZincLig',exptName
        
        lowef = LowEDir+ exptName + '_lowE.csv'
        if NBestLig == None:
            ligTbl, batch2ligTbl = loadBestLig(lowef)
        else:
            ligTbl, batch2ligTbl = loadBestLig(lowef,maxLig=NBestLig)
    
        print 'analFAAH_3: Exp=%s: NLig=%d' % (exptName,len(ligTbl))
       
        L2FPickleFile =   L2FDir + exptName + '_lig2Frag.pkl'
        if os.path.exists(L2FPickleFile):
            # print 'L2FPickle exists; using',L2FPickleFile
            l2fTbl = cPickle.load(open(L2FPickleFile,'rb'))
        else:
            print  'analFAAH: %s no L2FPickleFile found!? run bldLig2frag()' % (exptName)
            continue
                              
        bstart = exptTbl[expt]['bstart']
        bend = exptTbl[expt]['bend']
  
        print 'R2FPickleFile not found, building...',R2FPickleFile
        errf = R2FDir + exptName + 'r2f_err.txt'
        r2fTbl = bldRLIF2frag(exptName,ligTbl,l2fTbl,errf,LigPathTbl,bstart,bend) # RLIF -> frag -> [ (zincid,fragIdx,latomFull,fragInfo) ]
                             
        print 'R2FPickleFile dumping to',R2FPickleFile
        cPickle.dump(r2fTbl, open(R2FPickleFile,'wb'))
        r2fFile = R2FDir + exptName + '_rlif2Frag.csv'
        rptR2FTbl(r2fTbl,r2fFile)
  
        elapTime = datetime.datetime.now() - begTime
        print '</analFAAH_3 %s %s sec>' % (exptName,elapTime.seconds)
       
#     sys.exit( 'bldRLIF2frag done.' )
    
    ############################################################
    
    print 'analFAAH: *PHASE 4: capture ligand atom coordinates associated w/ RLIF'
    
    if not os.path.isdir(LigCoordDir):
        print 'analFAAH_4: creating LigCoord directory',LigCoordDir
        os.makedirs(LigCoordDir)
     
    for expt in allExpt:
             
        exptName = bldExptStr(expt)
    
        LigCoordPickleFile =   LigCoordDir + exptName + '_ligCoord.pkl'
        if os.path.exists(LigCoordPickleFile):
            print 'analFAAH_4: ligCoord.pkl already exists',exptName
            continue
            
        ## NonZincLigTbl is experiment-specific
        ## needs to be kept for reloading along with pickeled
        config.NonZincLigTbl = {} # ligand name -> nonzincID 
        config.NZIdx2LigTbl = {} # nonzincID -> ligand name
        config.NNonZincLig = 0
        if config.RunName.startswith('DUDE'):
            nonzf = LowEDir + '%s_nonZincLig.csv' % (exptName)
        else:
            nonzf = LowEDir + 'nonZincLig.csv'
        if os.path.exists(nonzf):
            # cf. analFAAH_expt()
            reader = csv.DictReader(open(nonzf))
            for i,entry in enumerate(reader):
                # Ligand,NZIdx
                config.NonZincLigTbl[ entry['Ligand'] ] = int(entry['NZIdx'])
                config.NZIdx2LigTbl [ int(entry['NZIdx']) ] =  entry['Ligand']
            config.NNonZincLig = len(config.NonZincLigTbl)
            print 'analFAAH_4: %d NonZincLig loaded from %s' % (config.NNonZincLig,nonzf)
        else:
            print 'analFAAH_4: no NonZincLig',exptName
      
        lowef = LowEDir+ exptName + '_lowE.csv'
        if NBestLig == None:
            ligTbl, batch2ligTbl = loadBestLig(lowef)
        else:
            ligTbl, batch2ligTbl = loadBestLig(lowef,maxLig=NBestLig)
    
        R2FPickleFile = R2FDir + exptName + '_rlif2frag.pkl'
        r2fTbl = cPickle.load(open(R2FPickleFile,'rb')) # RLIF -> frag -> [ (zincid,fragIdx) ]
     
        ligCoordTbl = bldLigCoordFromPDBQT(exptName,r2fTbl,ligTbl)
        cPickle.dump(ligCoordTbl, open(LigCoordPickleFile,'wb'))
     
#     sys.exit( ('Done dumping ligCoord' ) )
    
    
#     ############################################################
     
    print 'analFAAH: *PHASE 5: bldR2FCtr_v2'
         
    if not os.path.isdir(R2FCDir):
        print 'analFAAH_5: creating R2FCDir directory',R2FCDir
        os.makedirs(R2FCDir)
         
#     if not os.path.isdir(PlotDir):
#         print 'analFAAH_5: creating PlotDir directory',PlotDir
#         os.makedirs(PlotDir)
     
    for expt in allExpt:
      
        exptName = bldExptStr(expt)
             
        # NB: existence of BOTH fragClust, plot used to indicate phase already run
             
#         if config.RunName.startswith('DUDE') and exptName.find('w1c0') != -1:
#             print 'analFAAH_5: using w0c0 fragSimClust for w1c0 in 1_PR_x3KF0_prAS'
#             fragClustf = FragSimDir + '1_PR_x3KF0_prASw0c0_AS_DD_fragClust.csv'
#         else:
#             fragClustf = FragSimDir + exptName + '_fragClust.csv'
                 
        R2FCPickleFile = R2FCDir + exptName + '_r2fc.pkl'    
        if os.path.exists(R2FCPickleFile):
            print 'analFAAH Phase 5: fragClust already exist',exptName
            continue
   
        begTime = datetime.datetime.now()
        begTimeStr = begTime.strftime('%y%m%d_%H%M%S')
        print '<analFAAH_5  %s %s>' % (exptName,begTimeStr)
         
        R2FPickleFile = R2FDir + exptName + '_rlif2frag.pkl'
        r2fTbl = cPickle.load(open(R2FPickleFile,'rb')) # RLIF -> frag -> [ (zincid,fragIdx,latomFull,fragInfo) ]
     
        ## NonZincLigTbl is experiment-specific
        ## needs to be kept for reloading along with pickeled
        config.NonZincLigTbl = {} # ligand name -> nonzincID 
        config.NZIdx2LigTbl = {} # nonzincID -> ligand name
        config.NNonZincLig = 0
        if config.RunName.startswith('DUDE'):
            nonzf = LowEDir + '%s_nonZincLig.csv' % (exptName)
        else:
            nonzf = LowEDir + 'nonZincLig.csv'
        if os.path.exists(nonzf):
            # cf. analFAAH_expt()
            reader = csv.DictReader(open(nonzf))
            for i,entry in enumerate(reader):
                # Ligand,NZIdx
                config.NonZincLigTbl[ entry['Ligand'] ] = int(entry['NZIdx'])
                config.NZIdx2LigTbl [ int(entry['NZIdx']) ] =  entry['Ligand']
            config.NNonZincLig = len(config.NonZincLigTbl)
            print 'analFAAH_5: %d NonZincLig loaded from %s' % (config.NNonZincLig,nonzf)
        else:
            print 'analFAAH_5: no NonZincLig',exptName
             
        lowef = LowEDir+ exptName + '_lowE.csv'
        if NBestLig == None:
            ligTbl, batch2ligTbl = loadBestLig(lowef)
        else:
            ligTbl, batch2ligTbl = loadBestLig(lowef,maxLig=NBestLig)
     
        print 'analFAAH_5: Exp="%s": NLig=%d' % (exptName,len(ligTbl))
     
        r2fcfile = R2FCDir + exptName + '_r2fc.csv'
        r2fcTbl = bldR2FCtr_v2(exptName,r2fTbl,ligTbl,r2fcfile)
           
        cPickle.dump(r2fcTbl, open(R2FCPickleFile,'wb'))
   
        elapTime = datetime.datetime.now() - begTime
        print '</analFAAH_5 %s %s sec>' % (exptName,elapTime.seconds)
     
    # sys.exit( 'bldRLIF2fragCtr  done.' )
    
#     ############################################################
    
    timeRptFile = SummRptDir + 'fileTimes.csv'
    collectTimes(exptTbl,timeRptFile)
        
#     ############################################################
  
    if not(config.RunName.startswith('SAMPL4') or config.RunName.startswith('DUDE')):
        sys.exit( 'analFAAH: classification only defined over SAMPL4, DUDE.' )
  
    print 'analFAAH: *PHASE 6: build ARFF'

    if config.RunName.startswith('SAMPL4'):
        ## NonZincLigTbl is experiment-specific
        ## needs to be kept for reloading along with pickeled
        config.NonZincLigTbl = {} # ligand name -> nonzincID 
        config.NZIdx2LigTbl = {} # nonzincID -> ligand name
        config.NNonZincLig = 0
        nonzf = LowEDir + 'nonZincLig.csv' 
        # cf. analFAAH_expt()
        reader = csv.DictReader(open(nonzf))
        for i,entry in enumerate(reader):
            # Ligand,NZIdx
            config.NonZincLigTbl[ entry['Ligand'] ] = int(entry['NZIdx'])
            config.NZIdx2LigTbl [ int(entry['NZIdx']) ] =  entry['Ligand']
        config.NNonZincLig = len(config.NonZincLigTbl)
        print 'analFAAH_6: SAMPL4 %d NonZincLig loaded from %s' % (config.NNonZincLig,nonzf)
 
        SAMPLETrueTbl = loadSAMPLTrueTbl(SummRptDir+'ANSWERS/final_binders.txt')
        activeIdxSet = set(SAMPLETrueTbl.keys())
        print 'analFAAH_6: SAMPL4 NActive=%d' % (len(activeIdxSet))
          
    # NB: DUDE's actives loaded on per-experiment basis, below
   
    if not os.path.isdir(ArffDir):
        print 'analFAAH_6: creating ArffDir directory',ArffDir
        os.makedirs(ArffDir)
 
    if not os.path.isdir(PlotDir):
        print 'analFAAH_6: creating PlotDir directory',PlotDir
        os.makedirs(PlotDir)
     
    allExpt = exptTbl.keys()
    allExpt.sort()
     
    FragMinLigFreq = 2 # used in faahAnal.ligRFC2arff()
         
    for expt in allExpt:
        exptNo,prot,recept,site,lib = expt
        exptName = bldExptStr(expt)
  
        # NB: arrffFile written in ligRFC2arff()
        arrffFile = ArffDir+('%s.arff' % (exptName))
           
        if os.path.exists(arrffFile):
            print 'analFAAH_6: arrf already exist',exptName
            continue
                           
        if config.RunName.startswith('DUDE'):
             
            ## NonZincLigTbl is experiment-specific
            ## needs to be kept for reloading along with pickeled
            config.NonZincLigTbl = {} # ligand name -> nonzincID 
            config.NZIdx2LigTbl = {} # nonzincID -> ligand name
            config.NNonZincLig = 0
            nonzf = LowEDir + '%s_nonZincLig.csv' % (exptName)
            # cf. analFAAH_expt()
            if os.path.exists(nonzf):
                reader = csv.DictReader(open(nonzf))
                for i,entry in enumerate(reader):
                    # Ligand,NZIdx
                    config.NonZincLigTbl[ entry['Ligand'] ] = int(entry['NZIdx'])
                    config.NZIdx2LigTbl [ int(entry['NZIdx']) ] =  entry['Ligand']
                config.NNonZincLig = len(config.NonZincLigTbl)
                print 'analFAAH_6: %d NonZincLig loaded from %s' % (config.NNonZincLig,nonzf)
            else:
                print 'analFAAH_6: no NonZincLig',exptName
                 
            activeFile = ActiveDir + ActiveNameTbl[exptName]+'.lst'
            ins = open(activeFile)
            activeLigIdx = []
            for line in ins.readlines():
                lig = line.strip()
                activeLigIdx.append(normLigand(lig))
            ins.close()
            activeIdxSet = set(activeLigIdx)
            print 'analFAAH_6: DUDE %s NonZinc=%d NActive=%d' % (exptName,config.NNonZincLig,len(activeIdxSet))
  
        lowef = LowEDir+ exptName + '_lowE.csv'
        if NBestLig == None:
            ligTbl, batch2ligTbl = loadBestLig(lowef)
        else:
            ligTbl, batch2ligTbl = loadBestLig(lowef,maxLig=NBestLig)
           
        cummTrueMaxE, fracNeg = analTrueEnergy(ProcDir, PlotDir, exptName, ligTbl, activeIdxSet)
        print 'analFAAH_6: DUDE %s cummTrueMaxE=%f FracNeg=%f (%d)' % \
            (exptName,cummTrueMaxE,fracNeg,fracNeg*len(ligTbl))
          
        R2FCPickleFile = R2FCDir + exptName + '_r2fc.pkl'
        r2fcTbl = cPickle.load(open(R2FCPickleFile,'rb'))
        # NB: only use ethresh on biggest runs
        if len(ligTbl) > 20000:
            # don't write dense ARFF
            # ligRFC2arff(exptName,activeIdxSet,ethresh=cummTrueMaxE)
            ligRFC2SpArff(exptName,r2fcTbl,activeIdxSet,ethresh=cummTrueMaxE)
        else:
            # don't write dense ARFF
            # ligRFC2arff(exptName,activeIdxSet)
            ligRFC2SpArff(exptName,r2fcTbl,activeIdxSet)
            
    # sys.exit('done creating ARRF')
  
#     ############################################################
