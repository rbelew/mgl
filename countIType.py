'''
Created on Dec 14, 2014

@author: rik
'''

import argparse
import csv
import cPickle
import json
import re
import socket
import sys

### shared utilities
# 2do: shared with faahAnal
InterTypes = ('hba', 'hbd', 'mtl','ppi','tpi','vdw')

AARE = '([A-Z]+)([0-9]+)'
AAREPat = re.compile(AARE)

def bldExptStr(exptKey):
    s = '_'.join(exptKey)
    return s 

def ranges2list(rangeList):
    'invert bldRanges()'
    l = []
    for e in rangeList:
        if type(e)==type(1): # singleton int
            l.append(e)
        else:
            l += range(e[0],e[1]+1)
    return l

def bldFeaturePrefix(chain,raa):
    'C_999A: just chain + nicely formatted POSAAname'

    m = AAREPat.match(raa)
    try:
        (aaname,aapos) = m.groups()
    except:
        aapos = '999'
        aaname = raa

    iapos = int(aapos)
    # assert iapos < 1000, "iapos > 1000; '%03d' won't work"
    if iapos >= 1000:
        # print 'huh',raa, iapos
        # hack 141001!
        prefix = '%s_%s%s' % (chain,'999',aaname)
        
    aapos2 = '%03d' % iapos
    prefix = '%s_%s%s' % (chain,aapos2,aaname)
    
    return prefix

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

def bldExptTbl(inf):
    'return (exptNo,prot,recept,site,lib) -> exptData: bstart bend sys lib protein receptor site'
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

### end of shared utilities
def countAllIType(exptTbl,faahParentDir,runType,outf):
    allExpt = exptTbl.keys()
    allExpt.sort()
    outs = open(outf,'w')
    outs.write('Expt,TotLig,NX,N999')
    for it in InterTypes:
        outs.write(',%s_Lig,%s_Occ' % (it,it))
    outs.write('\n')
    outs.close()
      
    for exptKey in allExpt:
        exptNo,prot,recept,site,lib = exptKey
        faahDir = faahParentDir + 'Exp%s/' % (exptNo)

        exptData = exptTbl[exptKey]
        exptStr = bldExptStr(exptKey)
        
        batchList = ranges2list( [(exptData['bstart'],exptData['bend'])] )

        print 'coutAllIType: analyzing features %s %d batches' % \
            (exptKey,len(batchList))
        
        iDistTbl = countITypeDist(faahDir,runType,exptStr,batchList)
        totLig = iDistTbl['NRcd']
        nx = iDistTbl['NX']
        n999 = iDistTbl['N999']
        outs = open(outf,'a')
        outs.write('%s,%d,%d,%d' % (exptStr,totLig,nx,n999))
        for i in range(len(InterTypes)):                
            outs.write(',%d,%d' % (iDistTbl[i][0],iDistTbl[i][1]))
        outs.write('\n')
        outs.close()

def countITypeDist(faahDir,runType,exptName,bnoList):
    '''returns
            itypeTbl: itypeIdx: [nlig,nocc]
                        'NRcd': total number ligands
                        'NX' : number of bad X residue codes
                        'N999': number of bad residue positions
            distribution of interaction types across experiment
        first is number of ligands having ANY itype; second is total number interactions
        NB: assumes itypeIdx indices constant
    '''
          
    nrcd = 0
    nmissf = 0
    n999 = 0
    nx = 0
    itypeTbl = {}
    for i in range(len(InterTypes)):
        itypeTbl[i] = [0,0]  #nlig, totInteract
    
    for ib,bno in enumerate(bnoList):
        inf = faahDir+('inter/%s_inter_%07d.json' % (runType,bno))
        try:
            inStr = open(inf)
            allInter = json.load(inStr)
        except:
            nmissf += 1
            continue

        
        for il,interInfo in enumerate(allInter):
            nrcd += 1
            
            # cf crawl_ADV.rptData_ADV()
            # [ [Expt,BatchNo,Recept,Lig, [IType,[InterEnum] ] ] ]
            (expt,batchNo,recept,ligand, interList) = interInfo
            # NB: don't need to worry about ligand naming here

            for interTypeList in interList:
                itypeIdx, interEnum = interTypeList
                itypeTbl[itypeIdx][0] += 1
                itypeTbl[itypeIdx][1] += len(interEnum)
                itype = InterTypes[itypeIdx]
                for iinfo in interEnum:
                    (rchain,raa,ratom,liname) = getInterDetails(itype,iinfo)
                    prefix = bldFeaturePrefix(rchain,raa)
                    aas = prefix.split('_')[1]
                    if aas.startswith('999'):
                        n999 += 1
                    if aas.endswith('X'):
                        nx += 1

    itypeTbl['NRcd'] = nrcd    
    itypeTbl['NX'] = nx
    itypeTbl['N999'] = n999
        
    print 'countITypeDist: Expt=%s NRcd=%d NMissFile=%d NX=%d N999=%d' % (exptName,nrcd,nmissf,nx,n999)
    return itypeTbl

# arg string ala:
# 141214
# PrAS_115-120_141214 ADV  --verbose


if __name__ == '__main__': 

    parser = argparse.ArgumentParser(description='countAllIType arguments')
    parser.add_argument('runName',type=str,help='toplevel analysis topic name')
    parser.add_argument('runType',type=str,help='ADEngine runtype (AD, ADV)')
    parser.add_argument("--verbose",  action="store_true",help="increase output verbosity")
    
    args, unknown = parser.parse_known_args()

    HostName = socket.gethostname()
    if HostName == 'mgl0':
        print 'running on mgl0, good!'
        WCGDataDir = '/export/wcg/crawl/' 
        SummRptDir = '/export/wcg/anal/%s/'  % (args.runName)
    
    elif HostName == 'mgl3':
        print 'running on mgl3, slow(:'
        WCGDataDir = '/mgl/storage/wcg/crawl/' 
        SummRptDir = '/mgl/storage/wcg/anal/%s/'  % (args.runName)
    
    elif HostName.startswith('hancock'):
        print 'running local on hancock'
        WCGDataDir = '/Data/sharedData/coevol-HIV/WCG/processed/sampl/'
        SummRptDir = '/Data/sharedData/coevol-HIV/WCG/anal/%s/'  % (args.runName)
        
    else:
        print 
        sys.exit( ('unknown host %s' % (HostName)) )
    
    if args.verbose:
        print 'countAllIType: arguments'
        # NB: args is a Namespace object; 
        argsDict = vars(args)
        for k,v in argsDict.items():
            print '\t%s = %s' % (k,v)
    
    if len(unknown)>0:
        print 'countAllIType: huh?! Unkown arguments=', unknown
        assert False # can't do break or return here!

    exptTbl = bldExptTbl(SummRptDir+'exptBatch.csv')
            
    countAllIType(exptTbl,WCGDataDir,args.runType,SummRptDir+'itype.csv')