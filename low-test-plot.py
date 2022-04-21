##############################################################
#
# python low-test-plot.py --RunNumber 5872-000 --nev 3 --ReadFile True/False  --FileName scintE --Phase 1/2
#
##############################################################

import os, re, sys, glob, math, random, pickle
import numpy as np
import ROOT
import argparse
from scipy import stats
from scipy import interpolate 
from scipy.optimize import curve_fit
from collections import defaultdict
def platter_wf(wvs,numSamples,minSample):
    minChannel = 152
    apdSignals= {}
    for j in range(74):
        wv = wvs.GetWaveform(j)
        channel = wv.fChannel
        if channel < minChannel:
            continue
        apd = np.zeros(numSamples)
        for k in range(0,numSamples):
            np.put(apd, k, wv.At(k+minSample))
        apdSignals[channel] = apd
    platterSignals = np.zeros(2*numSamples)
    for i in range(numSamples):
        sum1,sum2 = 0,0
        for chn in apdSignals:
            if chn < 189:
                sum1 += apdSignals[chn][i]
            else:
                sum2 += apdSignals[chn][i]
        np.put(platterSignals,i,sum1)
        np.put(platterSignals,i+numSamples,sum2)
    return platterSignals
        
def read(runNumber,i,filename,phase):
    minChannel = 152
    minSample = 0
    numSamples = 2048
    root_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_root/phase%i/%s/'%(phase,filename)
    t = ROOT.TChain('tree','tree')
    print root_path+"temp_full_%s_gan_%s.root" % (runNumber,filename)
    t.Add(root_path+"temp_full_%s_gan_%s.root" % (runNumber,filename))
    ed = ROOT.EXOEventData()
    t.SetBranchAddress('EventBranch',ed)

    t2 = ROOT.TChain('tree','tree')
    print root_path+"temp_full_%s_real.root" % runNumber
    t2.Add(root_path+"temp_full_%s_real.root" % runNumber)
    ed2 = ROOT.EXOEventData()
    t2.SetBranchAddress('EventBranch',ed2)
    print "Same nEvents?", t.GetEntries() == t2.GetEntries()

    recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/%s/'%(phase,filename)
    tr = ROOT.TChain('tree','tree')
    print recon_path+"recon_%s_gan_%s.root" % (runNumber,filename)
    tr.Add(recon_path+"recon_%s_gan_%s.root" % (runNumber,filename))
    edr = ROOT.EXOEventData()
    tr.SetBranchAddress('EventBranch',edr)

    tr2 = ROOT.TChain('tree','tree')
    print recon_path+"recon_%s_real.root" % runNumber
    tr2.Add(recon_path+"recon_%s_real.root" % runNumber)
    edr2 = ROOT.EXOEventData()
    tr2.SetBranchAddress('EventBranch',edr2)

    t.GetEntry(i)
    t2.GetEntry(i)

    wvs = ed.GetWaveformData()
    wvs.Decompress()
    wvs2 = ed2.GetWaveformData()
    wvs2.Decompress()

    apdSignals,apdSignals2 = np.zeros(74*numSamples),np.zeros(74*numSamples)
    for j in range(74):
        wv = wvs.GetWaveform(j)
        wv2 = wvs2.GetWaveform(j)
        channel = wv.fChannel
        if channel < minChannel:
            continue
        BL = 0.0
        for k in range(2048/3):
            BL += wv.At(k)
        BL /= (2048/3)
        for k in range(0,numSamples):
            np.put(apdSignals, (channel-minChannel)*numSamples+k, wv.At(k+minSample)-BL)

        BL = 0.0
        for k in range(2048/3):
            BL += wv2.At(k)
        BL /= (2048/3)
        for k in range(0,numSamples):
            np.put(apdSignals2, (channel-minChannel)*numSamples+k, wv2.At(k+minSample)-BL)

    platterSignals = platter_wf(wvs,numSamples,minSample)
    platterSignals2 = platter_wf(wvs2,numSamples,minSample)
    tr.GetEntry(i)
    tr2.GetEntry(i)

    nsc = edr.GetNumScintillationClusters()
    nsc2 = edr2.GetNumScintillationClusters()
    ncc = edr.GetNumChargeClusters()
    ncc2 = edr2.GetNumChargeClusters()

    esc = 0
    for j in range(nsc):
        sc = edr.GetScintillationCluster(j)
        esc += sc.fRawEnergy
    esc2 = 0
    for j in range(nsc2):
        sc = edr2.GetScintillationCluster(j)
        esc2 += sc.fRawEnergy
    for j in range(ncc):
        cc = edr.GetChargeCluster(j)
        pos = np.array([cc.fX,cc.fY,cc.fZ])
    for j in range(ncc2):
        cc = edr2.GetChargeCluster(j)
        pos2 = np.array([cc.fX,cc.fY,cc.fZ])
    with open("./pkls/wf_samples/low-test-sample-%s-%i.pkl"%(filename,i),"wb") as file:
        pickle.dump(apdSignals,file)
        pickle.dump(apdSignals2,file)
        pickle.dump(esc ,file) # / 8078.9 * 2615
        pickle.dump(esc2 ,file) # / 7010.8 * 2615
        pickle.dump(platterSignals,file)
        pickle.dump(platterSignals2,file)
        pickle.dump(pos,file)
        pickle.dump(pos2,file)
def plot_sample(s1,s2,nev,e1,e2,pos1,pos2):
    plt.rcParams['figure.figsize'] = [24, 10]
    plt.subplot(1,2,1)
    plt.title(r"GAN. Event %i. $E_s$ = %.2f [A.U.] pos = (%.2f,%.2f,%.2f)" % (nev,e1,pos1[0],pos1[1],pos1[2]))
    for i in range(len(s1)):
        plt.plot(s1[i] + 20 * i)
    plt.xlabel(r"time $[\mu s]$",fontsize=20)
    plt.ylabel(r"ADC values + channel number $\times$ 20",fontsize=20)
    plt.ylim([-200,1600])
    plt.xlim([0,2048])
    plt.subplot(1,2,2)
    plt.title(r"Real. Event %i. $E_s$ = %.2f [A.U.] pos = (%.2f,%.2f,%.2f)" % (nev,e2,pos2[0],pos2[1],pos2[2]))
    for i in range(len(s2)):
        plt.plot(s2[i] + 20 * i)
    plt.xlabel(r"time $[\mu s]$",fontsize=20)
    plt.ylabel(r"ADC values + channel number $\times$ 20",fontsize=20)
    plt.ylim([-200,1600])
    plt.xlim([0,2048])
    plt.savefig("pics/wf_samples/waveforms_test_%i.png"%(nev),bbox_inches = 'tight')
    plt.show()
def plot_platter(p1,p2,nev,e1,e2,pos1,pos2):
    plt.rcParams['figure.figsize'] = [24, 10]
    wf1,wf2 = np.sum(p1,axis=0),np.sum(p2,axis=0)
    maxy = max(wf1.max(),wf2.max()) + 200
    miny = min(wf1.min(),wf2.min()) - 200
    plt.subplot(1,2,1)
    plt.title(r"GAN. Event %i. $E_s$ = %.2f [A.U.] pos = (%.2f,%.2f,%.2f)" % (nev,e1,pos1[0],pos1[1],pos1[2]))
    #for i in range(len(p1)):
    #    plt.plot(p1[i] + 20 * i)
    plt.plot(wf1)
    plt.xlabel(r"time $[\mu s]$",fontsize=20)
    plt.ylabel(r"ADC values + channel number $\times$ 20[a.u.]",fontsize=20)
    plt.ylim([miny,maxy])
    plt.xlim([0,2048])
    plt.subplot(1,2,2)
    plt.title(r"Real. Event %i. $E_s$ = %.2f [A.U.] pos = (%.2f,%.2f,%.2f)" % (nev,e2,pos2[0],pos2[1],pos2[2]))
    #for i in range(len(p2)):
    #    plt.plot(p2[i] + 20 * i)
    plt.plot(wf2)
    plt.xlabel(r"time $[\mu s]$")
    plt.ylabel(r"Plane number + [a.u.]")
    plt.ylim([miny,maxy])
    plt.xlim([0,2048])
    plt.savefig("pics/wf_samples/waveforms_platter_%i.png"%(nev),bbox_inches = 'tight')
    plt.show()
def plot_single_channel(s1,s2,nev,e1,e2,pos1,pos2,chn):
    plt.rcParams['figure.figsize'] = [24, 10]
    plt.subplot(1,2,1)
    plt.title(r"GAN. Event %i. channel %i" %(nev,chn + 152))
    plt.plot(s1[chn])
    plt.xlabel(r"time $[\mu s]$")
    plt.ylabel(r"Channel number + [a.u.]")
    plt.ylim([s1[chn].min(),s1[chn].max()])
    plt.xlim([0,2048])
    plt.subplot(1,2,2)
    plt.title(r"Real. Event %i. channel %i" % (nev,chn + 152))
    plt.plot(s2[chn])
    plt.xlabel(r"time $[\mu s]$")
    plt.ylabel(r"Channel number + [a.u.]")
    plt.ylim([s1[chn].min(),s1[chn].max()])
    plt.xlim([0,2048])
    plt.savefig("pics/wf_samples/waveforms_test_%i.png"%(nev),bbox_inches = 'tight')
    plt.show()
def plot_from_wf(s1,s2,nev,e1,e2,pos1,pos2):
    plt.title(r"Real. Event %i. $E_s$ = %.2f [A.U.] pos = (%.2f,%.2f,%.2f)" % (nev,e2,pos2[0],pos2[1],pos2[2]))
    for i in range(len(s2)):
        line = s2[i] + 20 * i
        plt.plot(np.array(range(850)),line[:850],color='0.8')
        plt.plot(np.array(range(850,1200)),line[850:1200],color='k')
        plt.plot(np.array(range(1200,2048)),line[1200:],color='0.8')
    plt.xlabel(r"time $[\mu s]$",fontsize=20)
    plt.ylabel(r"ADC values + channel number $\times$ 20",fontsize=20)
    plt.ylim([-200,1600])
    plt.xlim([0,2048])
    plt.savefig("pics/wf_samples/waveforms_sample_%i.png"%(nev),bbox_inches = 'tight')
    plt.show()
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--RunNumber', type=str, nargs=1, default=5872)
    parser.add_argument('--nev', type=int, nargs=1, default=0)
    parser.add_argument('--FileName', type=str, nargs=1)
    parser.add_argument('--ReadFile', type=bool,default=False)
    parser.add_argument('--Phase', type=int,default=1)
    args = parser.parse_args()
    RunNumber = args.RunNumber[0]
    readfile = args.ReadFile
    FileName = args.FileName[0]
    Phase = args.Phase
    nev = args.nev[0]
    if readfile:
        import matplotlib as mpl
        mpl.use('Agg')
        import matplotlib.pyplot as plt
        with open("./pkls/wf_samples/low-test-sample-%s-%i.pkl"%(FileName,nev),"rb") as file:
            s1 = pickle.load(file)
            s2 = pickle.load(file)
            e1 = pickle.load(file)
            e2 = pickle.load(file)
            p1 = pickle.load(file)
            p2 = pickle.load(file)
            pos1 = pickle.load(file)
            pos2 = pickle.load(file)
        file.close()
        s1,s2 = s1.reshape(74,2048),s2.reshape(74,2048)
        p1,p2 = p1.reshape(2,2048),p2.reshape(2,2048)
        print s1.shape,s2.shape,p1.shape,p2.shape
        plot_from_wf(s1,s2,nev,e1,e2,pos1,pos2)
        plot_sample(s1,s2,nev,e1,e2,pos1,pos2)
        #for i in [11,26,39,53,73]:
        #    plot_single_channel(s1,s2,nev,e1,e2,pos1,pos2,i)
        plot_platter(p1,p2,nev,e1,e2,pos1,pos2)
    else:
        import ROOT
        ROOT.gSystem.Load("libEXOUtilities")
        ROOT.gSystem.Load("libEXOFitting")
        read(RunNumber,nev,FileName,Phase)
