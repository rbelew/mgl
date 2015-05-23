# -*- coding: utf-8 -*-
"""
Created on Thur May 21 19:31:09 2015

@author: santiagodn, rik
"""

# rikPlot4dns.py

import os, math, sys
from operator import itemgetter
import csv

import matplotlib as mpl
import __main__
mpl.use('Agg')
#import mpl.pylab as pl
from matplotlib import pylab as pl

def myAuc(fpr, tpr,maxFpr=1.0,reorder=True):
    if reorder:
        fpr2=[]
        tpr2=[]
        tList=[]
        # !!! preserve pairings !!!
        # populate a list of (fpr,tpr)
        for i,x in enumerate(fpr):
            if x<=maxFpr:
                tList.append( (x,tpr[i]) )
        # sort the populated list
        tList.sort(key=itemgetter(0), reverse=False)
        
        # extract re-oredered (fpr,tpr)
        for t in tList:
            fpr2.append(t[0])
            tpr2.append(t[1])
    else:
        fpr2=fpr
        tpr2=tpr
    if len(fpr2) != len(tpr2):
        return -1
    sum=0
    previous_x=0
    height=0
    for i,x in enumerate(fpr2):
        if x<=maxFpr:
            sum+=height*(x-previous_x)
            height=tpr2[i]
            previous_x=x
    return sum#/maxFpr

def myRIE( ligandList, activesList, alpha ):
    # Trunchon & Bayly, JCIM 2007, 47, 488-508, eq. 34
    # Schrodinger, enrichment.py
    N=len(ligandList)
    n=len(activesList)
    Ra=float(n)/N
    activesRanks=[]

    for rank, ligand in enumerate(ligandList):
        if ligand in activesList:
            activesRanks.append(rank)
    wSum=sum([math.exp(-1*alpha*r_i/N) for r_i in activesRanks])
    rie=wSum/(Ra*(1.0-math.exp(-1*alpha))/(math.exp(alpha/N)-1.0))
    return rie

def myBEDROC( rie, alpha, totalLigands, totalActives ):
    # Trunchon & Bayly, JCIM 2007, 47, 488-508, eq. 36
    # Schrodinger, enrichment.py
    Ra=float(totalActives)/totalLigands
    alphaRa=alpha*Ra
    bedroc=None
    frac1=None
    frac2=None
    try:
        frac1=Ra*math.sinh(alpha/2.0)/(math.cosh(alpha/2.0)-math.cosh(alpha/2.0-Ra*alpha))
        frac2=1.0/(1.0-math.exp(alpha*(1.0-Ra)))
        bedroc=rie*frac1+frac2
    except Exception, e:
        print e
        print '>>>','alpha',alpha,'Ra',Ra,'totalActives',totalActives,'totalLigands',totalLigands
    return bedroc, alphaRa

def analyzeSlist(ligands): 
    
    nligands    = len(ligands)              # was dat_tot
    nactive     = len(Actives)              # was totalActives
    ndecoys     = nligands - nactive   # was totalDecoys
    prand       = float(nactive)/nligands

    #roc_dat=[]
    ef_dat=[]
    tpr=[]
    fpr=[]
    # *** initial values before "guessing"
    roc_tp=0
    roc_tn=nligands-nactive
    roc_fp=0
    roc_fn=nactive
    roc_fpr=-1
    roc_tpr=-1
    
    #foundActives=0 # Tunchon & Bayly, JCIM 2007, 47, 488-508 # debug function to make sure actives existed in list
    sumActivesRanks=0
    #2=len(knownInh)
    #max_tpr=-1
    # *** initial values before "guessing"

    rie_dat={}
    bedroc_dat={}
    #myRIE( ligandList, activesList, alpha )
    alphaList=[10.0, 20.0,30.0,40.0,50.0,100.0,500.0]
    for alpha in alphaList:
        nowRie=myRIE( ligands, Actives, alpha)
        rie_dat[alpha]=nowRie
        #myBEDROC( rie, alpha, totalLigands,  )
        nowBEDROC=myBEDROC(nowRie,alpha,nligands,nactive)
        bedroc_dat[alpha]=nowBEDROC # (bedroc, alphaRa)
    
    ligFound = 0
    decFound = 0
    pLig=[]
    pDec=[]

    for rank,lig in enumerate(ligands):
        #update_roc=0
        if lig in Actives:
            sumActivesRanks+=(rank+1)
            ligFound+=1
    # NOTE: knownNonZeroes <= nonZeroes <= total
    # http://en.wikipedia.org/wiki/Receiver_operating_characteristic 
    # ==> roc curve = fpr vs. tpr      
            roc_tp+=1
            roc_fn-=1
        else:
            decFound+=1
            roc_fp+=1
            roc_tn-=1
        pLigFound=100.0*ligFound/nactive
        pLig.append(pLigFound)
        pDecFound=100.0*decFound/ndecoys
        pDec.append(pDecFound)
        psamp=float(roc_tp)/(roc_tp+roc_fp)
        samp=float(roc_tp+roc_fp)/(nligands)
        if prand != 0:
            ef_now = psamp/prand
        else:
            ef_now = 0
        ef_dat.append( (samp, ef_now, lig, rank) )
        roc_n  = roc_tn+roc_fp
        roc_p  = roc_tp+roc_fn
        if roc_n > 0 and roc_p>0:
            roc_spc= float(roc_tn)/roc_n  # specificity or true negative rate

            roc_fpr= 1-roc_spc     # fall-out or false positive rate
            roc_tpr= float(roc_tp)/roc_p  # sensitivity or true positive rate
            tpr.append( roc_tpr )
            fpr.append( roc_fpr )

    auacFrac= float(sumActivesRanks)/(float(nactive)*nligands)
    auac=1-auacFrac

    statDict = {'pLig': pLig,
                'pDec':  pDec,
                'fpr':  fpr,
                'tpr': tpr,
                'ef_dat': ef_dat,
                'auac':  auac,
                'rie_dat':  rie_dat,
                'bedroc_dat': bedroc_dat}
    
    return statDict

#def plot_ef5(ligSet,qProtocol,Kset,cpeTol,efDir,  \
#             fprD, tprD, ef_datD, pLigD, pDecD, typeD, \
#             fprS1,tprS1,ef_datS1,pLigS1,pDecS1,typeS1, \
#             fprS2,tprS2,ef_datS2,pLigS2,pDecS2,typeS2, \
#             fprS3,tprS3,ef_datS3,pLigS3,pDecS3,typeS3, \
#             fprS4,tprS4,ef_datS4,pLigS4,pDecS4,typeS4, ): # an3.py
def plot_ef5_rocAlt(plotName,lblList,confusion_data ):
    # confusion_data = [ ( pLig, pDec, fpr,tpr,ef_dat,datType, auac, rie_dat, bedroc_dat ) ]
    
    efDir = SummRptDir+'ef/'
    if not os.path.isdir( efDir ):
        print 'plot_ef5_rocAlt: creating EF directory', efDir
        os.makedirs( efDir )
        
    ef_root = efDir+'ef_'+plotName
    ef_file_00 = ef_root+'_00.png'
    ef_file_10 = ef_root+'_10.png'
    
    rocAltDir=SummRptDir +'roc_alt/'
    if not os.path.isdir( rocAltDir ):
        print 'plot_ef5_rocAlt: creating ROCAlt directory', rocAltDir
        os.makedirs( rocAltDir )
    roc_file3 = rocAltDir+'rocAlt_'+plotName+'-pFound.png'
    
    # plot ef
    color_index=0
    pl.clf()
    for lbl in lblList:
        samp=[]
        ef=[]
        step=[]
        ef_dat = confusion_data[lbl]['ef_dat']

        for i,thing in enumerate(ef_dat):
            step.append( i )
            samp.append( thing[0] )
            ef.append(thing[1])
        ef_max=max(ef[1:])
        ef_max_pindex=ef.index(ef_max)
        x_max=samp[ef_max_pindex]
        p_max=100*x_max
        
        pl.plot(samp,  ef,  label=lbl, color=Colors[color_index % NColors])
        pl.axvline(x=x_max,  linestyle='dotted', color=Colors[color_index % NColors], label='%s  EF_%0.2f = %0.2f (max.)'%(lbl,  p_max,  ef_max))
        color_index += 1

    pl.axhline(y=1, linestyle='dotted',color='k',label='Random')
    pl.xlim([0.0, 1.0])
    pl.xlabel('Fraction of Dockings')
    pl.ylabel('Enrichment Factor')
    pl.title('EF: %s'%(plotName))
    pl.legend(loc="upper right",prop={'size':10})
    #pl.gca().set_aspect('equal', adjustable='box')
    pl.savefig(ef_file_00)

    # plot ef (top10%)
    color_index=0
    pl.clf()
        
    for lbl in lblList:

        samp=[]
        ef=[]
        step=[]
        ef_dat = confusion_data[lbl]['ef_dat']
        
        for i,thing in enumerate(ef_dat):
            step.append( i )
            samp.append( thing[0] )
            ef.append(thing[1])
        ef_max=max(ef[1:])
        ef_max_pindex=ef.index(ef_max)
        x_max=samp[ef_max_pindex]
        p_max=100*x_max
        
        pl.plot(samp,  ef, label=lbl, color=Colors[color_index % NColors], marker='o', linestyle='',linewidth='0', markersize=1.0, markeredgecolor='none')
        pl.axvline(x=x_max,  linestyle='dotted', color=Colors[color_index % NColors], label='%s EF_%0.2f = %0.2f (max.)'%(lbl,  p_max,  ef_max))
        color_index += 1

    pl.axhline(y=1, linestyle='dotted',color='k',label='Random')
    pl.xlim([0.0, 0.1])
    pl.xlabel('Fraction of Dockings')
    pl.ylabel('Enrichment Factor')
    pl.title('Top 10p - EF: %s'%(plotName))
    pl.legend(loc="upper right",prop={'size':10})
    #pl.gca().set_aspect('equal', adjustable='box')
    pl.savefig(ef_file_10)

    # Plot ??? roc_file3, pFound
    color_index=0
    pl.clf()

    for lbl in lblList:
        samp=[]
        ef=[]
        step=[]
        lDec=[]
        lRand=[]
        lRandX=[]

        #convert proportion to log values (decoys)
        previousLDec=-1.75
        nowL=-1.75
        pDec = confusion_data[lbl]['pDec']
        pLig = confusion_data[lbl]['pLig']
        for p in pDec:
            if p == 0:
                lDec.append( previousLDec )
            else:
                previousLDec=nowL
                nowL=math.log10(p)
                lDec.append( math.log10(p) )
        pl.plot(lDec, pLig, label=lbl,color=Colors[color_index % NColors])
        color_index += 1


    pl.xlim([-1, 2])
    pl.ylim([0, 100])
    pl.set_xscale='log'
    pl.xlabel('LOG_10(%) Decoys Found')
    pl.ylabel('% Ligands Found')
    pl.title('Alternate ROC Plot: %s'%(plotName))
    pl.legend(loc="upper left")
    #convert proportion to log values (for random line)
    for p in range(-100,200,1):
        lRand.append(math.pow(10,p/100.0))
        lRandX.append( p/100.0 )
    pl.plot(lRandX, lRand, linestyle='dotted',color='k',label='Random')
    #pl.gca().set_aspect('equal', adjustable='box')
    pl.savefig(roc_file3)

#def plot_roc5(rec,qFieldSource,qProtocol,K,cpeTol,rocDir,
#              fprD,tprD,ef_datD, # an5.py
#              fprS1,tprS1,ef_datS1,typeS1,
#              fprS2,tprS2,ef_datS2,typeS2,
#              fprS3,tprS3,ef_datS3,typeS3,
#              fprS4,tprS4,ef_datS4,typeS4):

def plot_roc5(plotName,lblList,confusion_data):
                  
    rocDir = SummRptDir+'roc/'
    if not os.path.isdir( rocDir ):
        print 'plot_roc5: creating ROC directory', rocDir
        os.makedirs( rocDir )

    roc_root = rocDir+'roc_'+plotName
    
    roc_auc ={}
    
    for rocLevel in ROCLevelList:
        roc_auc[rocLevel]={}
        
    for lbl in lblList:
        for rocLevel in ROCLevelList:
            cMaxFpr=float(rocLevel)/100
            #print rocLevel, cMaxFpr
            fpr = confusion_data[lbl]['fpr']
            tpr = confusion_data[lbl]['tpr']
            roc_auc[rocLevel][lbl] = myAuc(fpr, tpr, maxFpr=cMaxFpr,reorder=True)

    for rocLevel in ROCLevelList:
        if rocLevel=='100':
            roc_file = roc_root+'_00.png'
        elif rocLevel=='10':
            roc_file = roc_root+'_10.png'
        elif rocLevel=='1':
            roc_file = roc_root+'_01.png'
        cMaxFpr=float(rocLevel)/100
        pl.clf()
        color_index=0
        for lbl in lblList:
            fpr = confusion_data[lbl]['fpr']
            tpr = confusion_data[lbl]['tpr']
            pl.plot(fpr, tpr, label='%s (area = %0.4f)' % (lbl,roc_auc[rocLevel][lbl]), color=Colors[color_index % NColors])
            color_index += 1
    
        pl.plot([0, 1], [0, 1], 'k--')
        pl.xlim([0.0, cMaxFpr])
        pl.ylim([0.0, 1.0])

        pl.xlabel('Fraction of Decoys Found')
        pl.ylabel('Fraction of Actives Found')
        pl.title('ROC: %s (%s%%)'%(plotName, rocLevel))
        pl.legend(loc="upper left",prop={'size':10})
        #mpl.pylab.show()
        pl.savefig(roc_file)

    return roc_auc
    
def get_ligands(file_name): # get ligands as ranked by Rik's score
    ligand_name_field = 0
    ligands=[]
    fs=open(file_name, 'r')
    flines=fs.readlines()
    fs.close()
    for iligand,ligand in enumerate(flines):
        if iligand>0:
            ligands.append(ligand.split(',')[ligand_name_field])
    return ligands

def get_docking_rankings(file_name):
    ligand_name_field = 0
    docking_score_field = 5
    docking_data=[]
    dockings    =[]
    fs=open(file_name, 'r')
    flines=fs.readlines()
    fs.close()
    for iligand,ligand in enumerate(flines):
        if iligand>0:
            docking_data.append( (ligand.split(',')[ligand_name_field],float(ligand.split(',')[docking_score_field])) )
    docking_data.sort(key=itemgetter(1), reverse=False)
    for datum in docking_data:
        #print datum[1],
        dockings.append( datum[0] )
    #print ''
    return dockings
  
def getLigands2(inf): 
    reader = csv.DictReader(open(inf))        
    ligands=[]
    for i,entry in enumerate(reader):
        # Ligand,Actual,Predict,PrTrue,Err,E,FPRate,TPRate
        try:
            ligands.append(entry['Ligand'])
        except Exception, e:
            print 'huh',e 
            print
    return ligands

def getDNSRankings(inf):

    ligands = []
    allTbl = [ {} for i in range(4)]
    ligColNum=0
    colNum = [2,3,4,5]

    reader = csv.reader(open(inf))
    for row in reader:
        lig = row[ligColNum]
        ligands.append(lig)
        for i in range(4):
            allTbl[i][lig] = float(row[colNum[i]])
        
    allDataTbl = {}
    for i in range(4):
        lbl = 's%d' % (i+1)
        newList = ligands[:]
        newList.sort(key=lambda lig: allTbl[i][lig],reverse=True)
        allDataTbl[lbl] = newList
        
    return allDataTbl
    
def get_ipa_rankings(file_name):
    fs=open(file_name, 'r')
    flines=fs.readlines()
    fs.close()
    
    s1_data=[]
    s2_data=[]
    s3_data=[]
    s4_data=[]
    
    ligand_name_field=0
    s1_field=2
    s2_field=3
    s3_field=4
    s4_field=5
    #          0                   1    2      3    4      5       6        7
    #ZINC00037275prasD_out_Vina_VS,0,43.310,4.671,0.701,0.021,-8.900000,-0.405000

    for iligand,ligand in enumerate(flines):
        s1_data.append( (ligand.split(',')[ligand_name_field], float(ligand.split(',')[s1_field])) )
        s2_data.append( (ligand.split(',')[ligand_name_field], float(ligand.split(',')[s2_field])) )
        s3_data.append( (ligand.split(',')[ligand_name_field], float(ligand.split(',')[s3_field])) )
        s4_data.append( (ligand.split(',')[ligand_name_field], float(ligand.split(',')[s4_field])) )
    
    s1_data.sort(key=itemgetter(1), reverse=True)
    s2_data.sort(key=itemgetter(1), reverse=True)
    s3_data.sort(key=itemgetter(1), reverse=True)
    s4_data.sort(key=itemgetter(1), reverse=True)
    
    s1_list=[]
    s2_list=[]
    s3_list=[]
    s4_list=[]
    
    for datum in s1_data:
        #print datum[1],
        s1_list.append( datum[0] )
    for datum in s2_data:
        #print datum[1],
        s2_list.append( datum[0] )
    for datum in s3_data:
        #print datum[1],
        s3_list.append( datum[0] )
    for datum in s4_data:
        #print datum[1],
        s4_list.append( datum[0] )

    return s1_list, s2_list, s3_list, s4_list

if __main__:

    import socket
    HostName = socket.gethostname()
    if HostName == 'mgl3':
        print 'running on mgl3'
        BaseDir = '/mgl/storage/wcg/'
    
    elif HostName.startswith('hancock'):
        print 'running local on hancock'
        BaseDir = '/Data/sharedData/coevol-HIV/WCG/'
   
    else:
        print 
        sys.exit( ('unknown host %s' % (HostName)) )
        
    RunName = 'iniTst'

    SummRptDir = BaseDir + 'anal/%s/plots/'  % (RunName)
    DataDir = SummRptDir + 'dat/'
     
    # necessary actives list
    actives_file = DataDir + 'pras.lst'

    Actives = []
    # get actives
    fs=open(actives_file,'r')
    for active in fs.readlines():
        Actives.append(active.strip())
    fs.close()
    nactive=len(Actives)
    print 'rikPlot4dns: %d Actives read from %s' % (nactive,actives_file)

    # NB: original argv ordering of dataSrc and their labels maintained
    lblList = []
    dataSourceList = []
    allLigLists = {}
    confusion_data={}

    for dataFile in sys.argv[1:]:
        
        dataSrc = dataFile[:-4] # drop '.csv'
        # 2do-HACK: use "_dns" to flag DNS's result files
        if dataSrc.find('_dns') != -1:
            allEvalTbl = getDNSRankings(DataDir+dataFile) # elbl -> ligands
            # original DNS, for comparison
            s1_list, s2_list, s3_list, s4_list = get_ipa_rankings(DataDir+dataFile)
            allElbls = allEvalTbl.keys()
            allElbls.sort()
            for elbl in allElbls:
                if elbl in allLigLists:
                    print 'rikPlot4dns: duplicate data label1?!',elbl
                    continue
                
                dataSourceList.append(dataSrc)
                # NB: this breaks elbl into characters, then adds each character
                # lblList += elbl  
                lblList.append(elbl)
                
                ligands = allEvalTbl[elbl]
                allLigLists[elbl] = ligands
                statDict = analyzeSlist(ligands)
                confusion_data[elbl] = statDict
                 
        else:
            # 2do-HACK: use "_lbl.csv" to pass lbl
            rbpos = dataSrc.rfind('_')
            lbl = dataSrc[rbpos+1:]
            if lbl in allLigLists:
                print 'rikPlot4dns: duplicate data label2?!',lbl
                continue
            dataSourceList.append(dataSrc)
            lblList += lbl
            ligands = getLigands2(DataDir+dataFile)
            # original DNS, for comparison
            ligands_list = get_ligands(DataDir+dataFile)
            dockings_list = get_docking_rankings(DataDir+dataFile)

            allLigLists[lbl] = ligands
            statDict = analyzeSlist(ligands)
            confusion_data[lbl] = statDict
         
    ## Confirm that all ligand lists are coextensive
    allLig = set()
    for il,lbl in enumerate(lblList):
        allLig = allLig | set(allLigLists[lbl])
    ntotlig = len(allLig)
    print 'rikPlot4dns: NDataSets=%d NLigands=%d Labels=%s' % (len(lblList),ntotlig,lblList)
    if not all([len(allLigLists[lbl])==ntotlig for lbl in lblList]):
        print 'rikPlot4dns: Ligand lists mismatch?!  Union=%d' % (ntotlig)

        print 'Lbl\t% 25s\tNLig\tNMiss' % ('DataSource')
        for il,lbl in enumerate(lblList):
            if len(allLigLists[lbl]) != ntotlig:
                missingLig = allLig - set(allLigLists[lbl])
                print '%s\t%s\t%d\t%d' % \
                (lbl,dataSourceList[il],len(allLigLists[lbl]),len(missingLig))
                ## Extend lists with any missing ones
                
                allLigLists[lbl] += list(missingLig)
            
        # Try again
        if not all([len(allLigLists[lbl])==ntotlig for lbl in lblList]):   
            sys.exit( 'rikPlot4dns: Ligand lists STILL mismatch?!' )
    

    # output data
    ef_fractionList = [0.001, 0.002, 0.01, 0.05, 0.1, 0.2]
    ROCLevelList=['100','10','1'] 
    
    Colors=['r','k','b','m','g','c','o']
    NColors=len(Colors)
            
    # plot and calculate/process

    plot_ef5_rocAlt('initTst',lblList,confusion_data)

    roc_auc_list = plot_roc5('initTst',lblList,confusion_data)
    
#     datOutFile = lig_set+'_dat.csv'
#     datOutFileCount=0
#     datHeader='LigSet,Score,Data,SubData,Value\n'
#     
#     while os.path.isfile(datOutFile):
#         datOutFileCount+=1
#         datOutFile = lig_set+'_dat-'+str(datOutFileCount)+'.csv'
#     
#     ds=open(datOutFile,'w')
#     ds.write(datHeader)
#     
#     for confusion_key in confusion_data:
#         pLig, pDec, fpr,tpr,ef_dat,datType, auac, rie_dat, bedroc_dat = confusion_data[confusion_key]
#     
#         #AUAC
#         auac_line='%s,%s,%s,%s,%.4f\n' % (lig_set,datType,'auac','-',auac)
#         ds.write(auac_line)
#         
#         #AUC
#         for rocLevel in rocLevelList:
#             auc_line='%s,%s,%s,%s,%.4f\n' % (lig_set,datType,'auc',str(rocLevel),roc_auc_list[rocLevel][datType])
#             ds.write(auc_line)
#         
#         #BEDROC
#         for alpha in alphaList:
#             #print 'bedroc_dat[int(alpha)]',bedroc_dat[int(alpha)] = (bedroc,alphaRa)
#             bedroc_line='%s,%s,%s,%s,%.4f\n' % (lig_set,datType,'bedroc',str(alpha),bedroc_dat[int(alpha)][0])
#             ds.write(bedroc_line)
#     
#         # EF
#         # ef_dat = [ (samp, ef_now, thing, rank) ]
#         get_ef={}
#         for ef_fraction in ef_fractionList:
#             get_ef[ef_fraction]={}
#             get_ef[ef_fraction]['flag']=1
#         past_ef=0
#         now_ef=0
#         for ef_data in ef_dat:
#             samp, ef_now, thing, rank = ef_data
#             for ef_fraction in ef_fractionList:
#                 if get_ef[ef_fraction]['flag']:
#                     if samp > ef_fraction:
#                         #get_ef[ef_fraction]['value']=ef_now
#                         ef_line='%s,%s,%s,%s,%.4f\n' % (lig_set,datType,'ef',str(ef_fraction),ef_now)
#                         ds.write(ef_line)
#                         get_ef[ef_fraction]['flag'] =0
#         del get_ef
#         
#         #RIE
#         for alpha in alphaList:
#             rie_line='%s,%s,%s,%s,%.4f\n' % (lig_set,datType,'rie',str(alpha),rie_dat[int(alpha)])
#             ds.write(rie_line)
#     
#     ds.close()