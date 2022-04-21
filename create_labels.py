#################################################################################
#
# python create_labels.py --RunNumber 5872 S8 4769; S2 7050 wCo 5349 --FileNumber 0/1/2
#
# ncl=1 & nsc=1 & e_q (fPurityCorrectedEnergy)>600 & e_scint (fRawEnergy) > 0
# & radius<162 & abs(z)<180 & wf_data.fNumSamples=2048
#################################################################################

import csv
import ROOT
import numpy as np 
import math
import cPickle as pickle
import argparse
def save_real_wf(apdSignals,nEvents,runNumber,fileNumber):
    outDir = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_wf/real/phase2'
    option = 'a'
    if nEvents == 0:
        option = 'w'

    row = ''
    for k in range(len(apdSignals)):
        row += '%f'%apdSignals[k]
        if not k == len(apdSignals)-1:
            row += ','
    row += '\n'

    fname = outDir+'/real_waveform'
    f = open('%s_%i-%03i.csv'%(fname,runNumber,fileNumber),option)
    f.write(row)
    f.close()
def GetAPDSignals(waveformData, t_scint, subtractBaseline):
    numAPD = 74
    minChannel = 152
    minSample = 850
    numSamples = 350
    maxSample = minSample+numSamples
    nsamples = waveformData.fNumSamples
    nWF = waveformData.GetNumWaveforms()

    if t_scint > -180:
        minSample += int(t_scint)
        maxSample += int(t_scint)

    channelData = np.zeros((numAPD,numSamples))
    channels = np.zeros(numAPD)
    sumData = np.zeros(numSamples)

    apdSignals = np.zeros(numAPD*numSamples)

    if not nsamples == 2048:
        return apdSignals

    if t_scint < -180:
        return apdSignals

    if maxSample > 2048:
        return apdSignals

    n = 0
    for i in range(nWF):
        wf = waveformData.GetWaveform(i)
        ch = wf.fChannel

        if ch < minChannel:
            continue

        wf.Decompress()

        BL = 0.0
        for k in range(nsamples/3):
            BL += wf.At(k)
        BL /= (nsamples/3)

        for k in range(0,numSamples):
            if subtractBaseline:
                np.put(apdSignals, (ch-minChannel)*numSamples+k, wf.At(k+minSample)-BL)
            else:
                np.put(apdSignals, (ch-minChannel)*numSamples+k, wf.At(k+minSample))
            np.put(sumData, k, sumData[k]+wf.At(k+minSample)-BL)

        n += 1

    maxBin = np.argmax(sumData)
    apdMax = np.max(sumData)
    rms = 0
    for k in range(numSamples/3):
        rms += sumData[k]*sumData[k]
    rms /= numSamples/3
    rms = math.sqrt(rms)
    
    #if apdMax < 5*rms:
    #   return np.zeros(numAPD*numSamples)

    return apdSignals
def PullData(runNumber,fileNumber,real_wfs = False):
    print("*************************************************************")
    print("RunNumber:       %i"%runNumber)
    print("*************************************************************")
    if real_wfs == True:
        subtractBaseline = True
    
    waveformFiles = '/nfs/slac/g/exo_data[wildcard]/exo_data/data/WIPP/root/[runNumber]/'
    preprocessedFiles = "/nfs/slac/g/exo_data[wildcard]/exo_data/data/WIPP/masked/[runNumber]/"
    
    ################# Open waveform files ##############################
    print 'Opening waveform files in directory %s'%(waveformFiles.replace('[runNumber]','%i'%runNumber)) 
    tWF = ROOT.TChain('tree','tree')
    for i in range(2,8):
        wfFile = waveformFiles.replace('[wildcard]','%i'%i)+'run0000%i-%03i.root'%(runNumber,fileNumber) 
        print wfFile.replace('[runNumber]','%i'%runNumber)
        tWF.Add(wfFile.replace('[runNumber]','%i'%runNumber))
    
    raw_ED = ROOT.EXOEventData()
    tWF.SetBranchAddress('EventBranch',raw_ED)
    tWF.BuildIndex("EventBranch.fRunNumber","EventBranch.fEventNumber")
    print "tWF events:",tWF.GetEntries()
    ################# Open masked files ##############################
    print 'Opening masked files in directory %s'%(preprocessedFiles.replace('[runNumber]','%i'%runNumber))
    tProcessed = ROOT.TChain('tree','tree')
    for i in range(2,8):
        MFile = preprocessedFiles.replace('[wildcard]','%i'%i)+'masked0000%i-%03i.root'%(runNumber,fileNumber)
        print MFile.replace('[runNumber]','%i'%runNumber)
        tProcessed.Add(MFile.replace('[runNumber]','%i'%runNumber))

    rec_ED = ROOT.EXOEventData()
    tProcessed.SetBranchAddress('EventBranch',rec_ED)
    tProcessed.BuildIndex("EventBranch.fRunNumber","EventBranch.fEventNumber")
    print "tProcessed events:",tProcessed.GetEntries()
    '''
    waveformFiles = '/nfs/slac/g/exo_data4/exo_data/data/WIPP/root/5872/run00005872-000.root'
    preprocessedFiles = '/nfs/slac/g/exo_data3/exo_data/data/WIPP/masked/5872/masked00005872-000.root' 
    
    print('Opening preprocessed file %s'%(preprocessedFiles))
    fProcessed = ROOT.TFile(preprocessedFiles,'READ')
    tProcessed = fProcessed.Get('tree')#dataTree')
    print("tProcessed events:",tProcessed.GetEntries())

    print('Opening waveform files %s'%(waveformFiles))
    tWF = ROOT.TChain('tree','tree')
    tWF.Add(waveformFiles)
    print("tWF events:",tWF.GetEntries())

    rec_ED = ROOT.EXOEventData()
    tProcessed.SetBranchAddress('EventBranch',rec_ED)

    raw_ED = ROOT.EXOEventData()
    tWF.SetBranchAddress('EventBranch',raw_ED)
    '''
    n = tProcessed.GetEntries()
    ni = 0
    labels = []
    picked = []
    for i in range(n):
        #if ni > 3000: break
        tProcessed.GetEntry(i)
        tWF.GetEntry(i)
        wf_data = raw_ED.GetWaveformData()

        ncl = rec_ED.GetNumChargeClusters()
        eq = rec_ED.GetTotalPurityCorrectedEnergy()
        nsc = rec_ED.GetNumScintillationClusters()
        if i % 1000 == 0: print 'i=%i, ni=%i'%(i,ni)
        if ncl<1 or nsc!=1 or wf_data.fNumSamples!=2048:
            continue
        if ncl == 1:
            cc = rec_ED.GetChargeCluster(0)
            sc = rec_ED.GetScintillationCluster(0)
            x,y,z,e_charge,e_scint = cc.fX,cc.fY,cc.fZ,cc.fPurityCorrectedEnergy,sc.fRawEnergy
            t_scint = sc.fTime / 1000 - rec_ED.fEventHeader.fTriggerOffset
            if x ** 2 + y ** 2 < 162 ** 2 and abs(z) < 180 and e_charge > 600 and e_scint > 0:
                labels.append([x,y,z,e_charge,e_scint])
                picked.append(i)
                if real_wfs == True:
                    apdSignals = GetAPDSignals(wf_data, t_scint, subtractBaseline)
                    save_real_wf(apdSignals,ni,runNumber,fileNumber)
                ni += 1

    #print picked[:100]
    labels = np.array(labels)
    print "nEvents: %i" %labels.shape[0]
    label_path = "/nfs/slac/g/exo/shaolei/GAN-test/label/"
    with open(label_path + "labels_%s-%03i.pkl" % (runNumber,fileNumber),"wb") as file:
        pickle.dump(labels,file,protocol=1)
if __name__ == "__main__":
    ROOT.gSystem.Load("libEXOUtilities")
    parser = argparse.ArgumentParser()
    parser.add_argument('--RunNumber', type=int, nargs=1, default=5972)
    parser.add_argument('--FileNumber', type=int, nargs=1, default=0)

    args = parser.parse_args()
    RunNumber = args.RunNumber[0]
    FileNumber = args.FileNumber[0]
    PullData(RunNumber,FileNumber,True)
