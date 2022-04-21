##############################################################
#
# python Amp_Z.py --FileName unflatten --ReadFile True/False --Phase 1/2
#
##############################################################

import os, re, sys, glob, math, random, pickle, csv
import numpy as np
import ROOT
import argparse
import matplotlib as mpl 
#mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
from scipy import interpolate 
from scipy.optimize import curve_fit
from collections import defaultdict
import ReconEvents as RE
def xy_to_phi(x,y):
    r = np.sqrt(x**2+y**2)
    if y >= 0:
        phi = np.arccos(x/r)
    else:
        phi = -np.arccos(x/r)
    return phi
def read(runNumbers,filename,phase):
    recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/%s'%(phase,filename) # gan_lightmap_recon/ recon_with_lightmap/lightmap_2017_0nu-base64_v2
    #recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/recon_with_lightmap/'%(phase)
    t = ROOT.TChain('tree','tree')
    for runNumber in runNumbers:
        t.Add(recon_path+"/recon_%s*_gan_%s.root" % (runNumber,filename))
        #print 'Add file:',runNumber
    ed = ROOT.EXOEventData()
    t.SetBranchAddress('EventBranch',ed)

    t2 = ROOT.TChain('tree','tree')
    for runNumber in runNumbers:
        t2.Add(recon_path+"/recon_%s*_real.root" % runNumber)
        #print 'Add file:',runNumber
    ed2 = ROOT.EXOEventData()
    t2.SetBranchAddress('EventBranch',ed2)
    print "Same nEvents?", t.GetEntries() == t2.GetEntries(),'gan:',t.GetEntries(),'real:',t2.GetEntries()

    wrong_event = []
    wrong_event_2 = []
    z,z2 = [],[]
    r,r2 = [],[]
    phi,phi2 = [],[]
    gan_amps,real_amps = [],[]
    gan_se,real_se = [],[]
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        t2.GetEntry(i)
        nev = ed.fEventNumber
        nev2 = ed2.fEventNumber

        nsc = ed.GetNumScintillationClusters()
        nsc2 = ed2.GetNumScintillationClusters()
        if nsc != nsc2:
            wrong_event.append([i,nev,nev2,nsc,nsc2])
        if nsc == 1:
            sc = ed.GetScintillationCluster(0)
            napd = sc.GetNumAPDSignals()
            ncc = sc.GetNumChargeClusters()
            if ncc == 1:
                cc = sc.GetChargeClusterAt(0)
                try:
                    cc.fPurityCorrectedEnergy
                except:
                    #print i
                    continue
                if cc.fPurityCorrectedEnergy > 2615+200 or cc.fPurityCorrectedEnergy < 2615-200: continue
                '''
                #wv1 = sc.GetAPDSignal(ROOT.EXOAPDSignal.kPlaneFit,1)
                #wv2 = sc.GetAPDSignal(ROOT.EXOAPDSignal.kPlaneFit,2)
                #GetCountsSumOnAPDPlane(EXOMiscUtil::kNorth)
                #wv1 = sc.GetPlaneOneSignal()
                #wv2 = sc.GetPlaneTwoSignal()
                plane1,plane2 = 0,0
                if wv1: plane1 = wv1.fRawCounts
                if wv2: plane2 = wv2.fRawCounts
                '''
                gan_se.append(sc.fRawEnergy)
                plane1,plane2 = sc.GetCountsSumOnAPDPlane(ROOT.EXOMiscUtil.kNorth),sc.GetCountsSumOnAPDPlane(ROOT.EXOMiscUtil.kSouth)
                amps = [plane1,plane2]
                gan_amps.append(amps)
                z.append(cc.fZ)
                r.append(cc.fX ** 2 + cc.fY ** 2)
                phi.append(xy_to_phi(cc.fX,cc.fY))
        if nsc2 == 1:
            sc = ed2.GetScintillationCluster(0)
            napd = sc.GetNumAPDSignals()
            ncc2 = sc.GetNumChargeClusters()
            if ncc2 == 1:
                cc = sc.GetChargeClusterAt(0)
                try:
                    cc.fPurityCorrectedEnergy
                except:
                    continue
                if cc.fPurityCorrectedEnergy > 2615+200 or cc.fPurityCorrectedEnergy < 2615-200: continue
                '''
                wv1 = sc.GetPlaneOneSignal()
                wv2 = sc.GetPlaneTwoSignal()
                plane1,plane2 = 0,0
                if wv1: plane1 = wv1.fRawCounts
                if wv2: plane2 = wv2.fRawCounts
                '''
                real_se.append(sc.fRawEnergy)
                plane1,plane2 = sc.GetCountsSumOnAPDPlane(ROOT.EXOMiscUtil.kNorth),sc.GetCountsSumOnAPDPlane(ROOT.EXOMiscUtil.kSouth)
                amps = [plane1,plane2]
                real_amps.append(amps)
                z2.append(cc.fZ)
                r2.append(cc.fX ** 2 + cc.fY ** 2)
                phi2.append(xy_to_phi(cc.fX,cc.fY))
                
    with open("./pkls/gan_real_Amp_Z.pkl","wb") as file:
        pickle.dump(np.array(gan_amps),file)
        pickle.dump(np.array(real_amps),file)
        pickle.dump(np.array(z),file)
        pickle.dump(np.array(z2),file)
        pickle.dump(np.array(r),file)
        pickle.dump(np.array(r2),file)
        pickle.dump(np.array(phi),file)
        pickle.dump(np.array(phi2),file)
        pickle.dump(np.array(gan_se),file)
        pickle.dump(np.array(real_se),file)
    print len(gan_amps),len(z)
    print len(real_amps),len(z2)
def get_mean_amp(amps,index,n):
    sum_amp = [[] for i in range(n)]
    for i,amp in zip(index,amps):
        try:
            sum_amp[i].append(amp)
        except:
            continue
    mean_amp,err_amp,std_amp = np.zeros(n),np.zeros(n),np.zeros(n)
    for i in range(n):
        if sum_amp[i]:
            np.put(mean_amp,i,np.mean(sum_amp[i]))
            np.put(err_amp,i,np.std(sum_amp[i]) / np.sqrt(len(sum_amp[i])))
            np.put(std_amp,i,np.std(sum_amp[i]))
    #mean_amp = [np.mean(i) for i in sum_amp]
    #err_amp = [np.std(i) / np.sqrt(len(i))for i in sum_amp]
    return mean_amp,err_amp,std_amp
def get_centers(bins,step):
    centers_list = []
    n = len(bins)-1
    for i in range(0,n):
        e_i = bins[i]+step/2.
        centers_list.append(e_i)
    return np.asarray(centers_list)
def plot_apd_phi(gan_amps,real_amps,gan_phi,real_phi,chn):
    left = -np.pi
    right = np.pi
    nbins = 8
    mybins = np.linspace(left,right,nbins+1)
    e_step = (right-left)/nbins
    if chn == 2:
        name = 'sum'
        outname = name
        gan_chn = [np.sum(i) for i in gan_amps]
        real_chn = [np.sum(i) for i in real_amps]
    else:
        name = str(chn)
        outname = str(chn+1)
        gan_chn = [i[chn] for i in gan_amps]
        real_chn = [i[chn] for i in real_amps]

    nGAN, binsGAN = np.histogram(gan_phi,bins=mybins)
    nReal,binsReal = np.histogram(real_phi,bins=mybins)
    index_gan = np.digitize(gan_phi,binsGAN) - 1
    index_real = np.digitize(real_phi,binsReal) - 1
    mean_amp_gan,err_amp_gan,std_amp_gan = get_mean_amp(gan_chn,index_gan,len(nGAN))
    mean_amp_real,err_amp_real,std_amp_real = get_mean_amp(real_chn,index_real,len(nReal))
    cent  = get_centers(mybins,e_step)

    plt.scatter(gan_phi,gan_chn,color='r',label = "GAN plane %s" % name)
    plt.scatter(real_phi,real_chn,color='b',label = "Real plane %s" % name)
    plt.title("fAPDs.fRawCounts")
    plt.xlabel(r"$\phi$ [rad]",fontsize=20)
    plt.ylabel(r"Amp Units",fontsize=20)
    plt.savefig("pics/phi_amp_gan_real_%s.png"%(outname),bbox_inches = 'tight')
    plt.show()
    plt.clf()

    #plt.errorbar(cent,mean_amp_real,yerr=err_amp_real,ecolor='lightgrey',capsize=0,elinewidth=2,fmt='o',color='b',markersize=4.,label=r'Real Amp plane %i' % chn)
    #plt.errorbar(cent,mean_amp_gan,yerr=err_amp_gan,ecolor='lightgrey',capsize=0,elinewidth=2,fmt='o',color='r',markersize=4.,label=r'GAN Amp plane %i' % chn)
    plt.plot(cent,mean_amp_gan,color='r',label = "GAN plane %s" % name)
    plt.fill_between(cent,mean_amp_gan-std_amp_gan,mean_amp_gan+std_amp_gan,color='r', alpha=0.2)
    plt.plot(cent,mean_amp_real,color='b',label = "Real plane %s" % name)
    plt.fill_between(cent,mean_amp_real-std_amp_real,mean_amp_real+std_amp_real, alpha=0.2)
    plt.legend()
    #Accuracy of axis labels tells people if you are a professional or a hack!
    plt.xlabel(r'$\phi$ [rad]',fontsize=20)
    plt.ylabel(r'Amp units',fontsize=20)
    plt.legend()
    
    #Set axes limits
    plt.axis((left,right,0,np.max(mean_amp_real)*1.75))
    plt.savefig("pics/phi_amp_map_%s.png"%(outname),bbox_inches = 'tight')
    plt.show()
    plt.clf()
def plot_apd_z(gan_amps,real_amps,gan_z,real_z,chn):
    left = -200
    right = 200
    nbins = 10
    mybins = np.linspace(left,right,nbins+1)
    e_step = (right-left)/nbins
    if chn == 2:
        name = 'sum'
        outname = name
        gan_chn = [np.sum(i) for i in gan_amps]
        real_chn = [np.sum(i) for i in real_amps]
    else:
        name = str(chn)
        outname = str(chn+1)
        gan_chn = [i[chn] for i in gan_amps]
        real_chn = [i[chn] for i in real_amps]

    nGAN, binsGAN = np.histogram(gan_z,bins=mybins)
    nReal,binsReal = np.histogram(real_z,bins=mybins)
    index_gan = np.digitize(gan_z,binsGAN) - 1
    index_real = np.digitize(real_z,binsReal) - 1
    mean_amp_gan,err_amp_gan,std_amp_gan = get_mean_amp(gan_chn,index_gan,len(nGAN))
    mean_amp_real,err_amp_real,std_amp_real = get_mean_amp(real_chn,index_real,len(nReal))
    cent  = get_centers(mybins,e_step)

    plt.scatter(gan_z,gan_chn,color='r',label = "GAN plane %s" % name)
    plt.scatter(real_z,real_chn,color='b',label = "Real plane %s" % name)
    plt.legend()
    plt.title("fAPDs.fRawCounts")
    plt.xlabel(r"Z [cm]",fontsize=20)
    plt.ylabel(r"Amp Units",fontsize=20)
    plt.savefig("pics/z_amp_gan_real_%s.png"%(outname),bbox_inches = 'tight')
    plt.show()
    plt.clf()

    #plt.errorbar(cent,mean_amp_real,yerr=err_amp_real,ecolor='lightgrey',capsize=0,elinewidth=2,fmt='o',color='b',markersize=4.,label=r'Real Amp plane %i' % chn)
    #plt.errorbar(cent,mean_amp_gan,yerr=err_amp_gan,ecolor='lightgrey',capsize=0,elinewidth=2,fmt='o',color='r',markersize=4.,label=r'GAN Amp plane %i' % chn)
    plt.plot(cent,mean_amp_gan,color='r',label = "GAN plane %s" % name)
    plt.fill_between(cent,mean_amp_gan-std_amp_gan,mean_amp_gan+std_amp_gan,color='r', alpha=0.2)
    plt.plot(cent,mean_amp_real,color='b',label = "Real plane %s" % name)
    plt.fill_between(cent,mean_amp_real-std_amp_real,mean_amp_real+std_amp_real, alpha=0.2)
    plt.legend()
    #Accuracy of axis labels tells people if you are a professional or a hack!
    plt.xlabel(r'Z [cm]',fontsize=20)
    plt.ylabel(r'Amp units',fontsize=20)
    plt.legend()
    
    #Set axes limits
    plt.axis((left,right,0,np.max(mean_amp_real)*1.75))
    plt.savefig("pics/z_amp_map_%s.png"%(outname),bbox_inches = 'tight')
    plt.show()
    plt.clf()
def plot_apd_r(gan_amps,real_amps,gan_r,real_r,chn):
    left = 0
    right = 180 ** 2
    nbins = 10
    mybins = np.linspace(left,right,nbins+1)
    e_step = (right-left)/nbins
    if chn == 2:
        name = 'sum'
        outname = name
        gan_chn = [np.sum(i) for i in gan_amps]
        real_chn = [np.sum(i) for i in real_amps]
    else:
        name = str(chn)
        outname = str(chn+1)
        gan_chn = [i[chn] for i in gan_amps]
        real_chn = [i[chn] for i in real_amps]

    nGAN, binsGAN = np.histogram(gan_r,bins=mybins)
    nReal,binsReal = np.histogram(real_r,bins=mybins)
    index_gan = np.digitize(gan_r,binsGAN) - 1
    index_real = np.digitize(real_r,binsReal) - 1
    mean_amp_gan,err_amp_gan,std_amp_gan = get_mean_amp(gan_chn,index_gan,len(nGAN))
    mean_amp_real,err_amp_real,std_amp_real = get_mean_amp(real_chn,index_real,len(nReal))
    cent  = get_centers(mybins,e_step)

    plt.scatter(gan_r,gan_chn,color='r',label = "GAN plane %s" % name)
    plt.scatter(real_r,real_chn,color='b',label = "Real plane %s" % name)
    plt.legend()
    plt.xlim([0,40000])
    plt.title("fAPDs.fRawCounts")
    plt.xlabel(r"$R^2$ $[cm^2]$",fontsize=20)
    plt.ylabel(r"Amp Units",fontsize=20)
    plt.savefig("pics/r_amp_gan_real_%s.png"%(outname),bbox_inches = 'tight')
    plt.show()
    plt.clf()

    #plt.errorbar(cent,mean_amp_real,yerr=err_amp_real,ecolor='lightgrey',capsize=0,elinewidth=2,fmt='o',color='b',markersize=4.,label=r'Real Amp plane %i' % chn)
    #plt.errorbar(cent,mean_amp_gan,yerr=err_amp_gan,ecolor='lightgrey',capsize=0,elinewidth=2,fmt='o',color='r',markersize=4.,label=r'GAN Amp plane %i' % chn)
    plt.plot(cent,mean_amp_gan,color='r',label = "GAN plane %s" % name)
    plt.fill_between(cent,mean_amp_gan-std_amp_gan,mean_amp_gan+std_amp_gan,color='r', alpha=0.2)
    plt.plot(cent,mean_amp_real,color='b',label = "Real plane %s" % name)
    plt.fill_between(cent,mean_amp_real-std_amp_real,mean_amp_real+std_amp_real, alpha=0.2)
    #Accuracy of axis labels tells people if you are a professional or a hack!
    plt.xlabel(r'$R^2$ $[cm^2]$',fontsize=20)
    plt.ylabel(r'Amp units',fontsize=20)
    plt.legend()

    #Set axes limits
    plt.axis((left,right,0,np.max(mean_amp_real)*1.75))
    plt.savefig("pics/r_amp_map_%s.png"%(outname),bbox_inches = 'tight')
    plt.show()
    plt.clf()
def plot_se_amp(se,amp,name=''):
    if len(se) != len(amp):
        print('Size not match.')
        return 0
    else:
        amp = np.sum(amp,axis=1)
        plt.scatter(se,amp,color='r',label = "%s Amp vs. Light E" % name)
        plt.legend()
        plt.xlim([2000,11000])
        plt.ylim([2000,11000])
        plt.xlabel('Amp Units')
        plt.ylabel('Light E (a.u.)')
        plt.savefig('pics/%s_amp_se.png'%name,bbox_inches='tight')
        plt.show()
def plot_amp_minus_energy(gan_se,real_se,gan_amp,real_amp):
    gan_sigma = gan_se-np.sum(gan_amps,axis=1)
    real_sigma = real_se-np.sum(real_amps,axis=1)
    mybins=np.linspace(gan_sigma.min(),gan_sigma.max(),40)
    plt.hist(gan_sigma,bins=mybins,label='GAN',alpha=0.5)
    plt.hist(real_sigma,bins=mybins,label='Real',alpha=0.5)
    plt.xlabel('fRawEnergy - fRawCounts')
    plt.legend()
    plt.show()
    outliers = [i for i in range(len(gan_sigma)) if np.abs(gan_sigma[i]) > 20]
    return outliers
def fit_gaussian(gans,reals,first,bounds,gan_first,gan_bounds,show=False,name = '',unit='keV'):
    #set the figure size
    plt.rcParams['figure.figsize'] = [12, 10]

    #create binning
    left = 500.
    right = bounds[1] + 200#12000
    nbins = 200#(right - left) / 20
    mybins = np.linspace(left,right,nbins+1)
    e_step = (right-left)/nbins

    #get bin content
    nreals, bins  = np.histogram(reals, bins=mybins)
    ngens, bins   = np.histogram(gans,  bins=mybins)
    
    err_r = RE.get_errs(nreals)
    err_g = RE.get_errs(ngens)
    cent  = RE.get_centers(mybins,e_step)

    #now actually plot
    #plt.step(cent,nlabels,'k',where='mid',linewidth=2) #try adding this for coarse binning - may look nicer
    plt.errorbar(cent,nreals,yerr=err_r,ecolor='lightgrey',capsize=0,elinewidth=2,fmt='o',color='b',markersize=4.,label=r'Light E (Reco on Real Waveforms)')
    plt.errorbar(cent,ngens,yerr=err_g,ecolor='lightgrey',capsize=0,elinewidth=2,fmt='o',color='r',markersize=4.,label=r'Light E (Reco on GAN Waveforms)')

    #Accuracy of axis labels tells people if you are a professional or a hack!
    plt.xlabel(r'E (%s)'%unit,fontsize=20)
    plt.ylabel(r'Counts/(%.0f %s)'%(e_step,unit),fontsize=20)

    #Set axes limits
    plt.axis((left,right,0,np.max(nreals)*1.75))

    #fit range
    left_r  = bounds[0]
    right_r = bounds[1]
    left_g  = gan_bounds[0]
    right_g = gan_bounds[1]
    cent_range_r, nrange_r, err_range_r = RE.get_range(cent,nreals,err_r,left_r,right_r)
    cent_range_g, nrange_g, err_range_g = RE.get_range(cent,ngens,err_g,left_g,right_g)
    xr = np.linspace(left_r,right_r,500)
    xg = np.linspace(left_g,right_g,500)

    #now fit and plot results
    popt_r, pcov = curve_fit(RE.gaussian, cent_range_r, nrange_r,first, maxfev=10000,sigma=err_range_r)
    perr_r = np.sqrt(np.diag(pcov))
    mychi2,ndof = RE.my_chisq(RE.gaussian,popt_r,cent_range_r,nrange_r,err_range_r)
    plt.plot(xr, RE.gaussian(xr, *popt_r), 'b-',linewidth=5.0,label=r'$\mu$=%5.1f$\pm$%3.1f %s, $\sigma$/$\mu$=%1.1f$\pm$%1.1f%%,  $\chi$$^2$/ndf=%.1f' % (popt_r[1],perr_r[1],unit,popt_r[2]/popt_r[1]*100.,perr_r[2]/popt_r[1]*100.,mychi2/ndof) )
    res_r = popt_r[2]/popt_r[1]*100.
    mu_r = popt_r[1]

    popt_g, pcov = curve_fit(RE.gaussian, cent_range_g, nrange_g,gan_first, maxfev=10000,sigma=err_range_g)
    perr_g = np.sqrt(np.diag(pcov))
    mychi2,ndof = RE.my_chisq(RE.gaussian,popt_g,cent_range_g,nrange_g,err_range_g)
    plt.plot(xg, RE.gaussian(xg, *popt_g), 'r-',linewidth=5.0,label=r'$\mu$=%5.1f$\pm$%3.1f %s, $\sigma$/$\mu$=%1.1f$\pm$%1.1f%%,  $\chi$$^2$/ndf=%.1f' % (popt_g[1],perr_g[1],unit,popt_g[2]/popt_g[1]*100.,perr_g[2]/popt_g[1]*100.,mychi2/ndof) )
    res_g = popt_g[2]/popt_g[1]*100.
    mu_g = popt_g[1]

    plt.legend(loc='upper left',numpoints = 1,frameon=True,prop={'size': 16})
    plt.gca().tick_params(labelsize=18)
    #plt.title('light E spectrum for %s'%name)

    plt.savefig('./pics/Th_Peak_%s.png'%name,bbox_inches = 'tight')
    if show:
        plt.show()
    else:
        plt.clf()
    return res_r,res_g,mu_r,mu_g

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #parser.add_argument('--RunNumber', type=int, nargs='+', default=5872)
    parser.add_argument('--FileName', type=str, nargs=1, default='unflatten')
    parser.add_argument('--ReadFile', type=bool,default=False)
    parser.add_argument('--Phase',type=int,default=1)
    args = parser.parse_args()
    #RunNumber = args.RunNumber
    readfile = args.ReadFile
    FileName = args.FileName[0]
    Phase = args.Phase
    print readfile
    RunNumber = []
    with open('phase%i_runset.csv'%Phase,'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if row[1] == 'wTh' and row[3] == 'test' and row[0] > 3717:
                RunNumber.append(row[0])
    if readfile:
        with open("./pkls/gan_real_Amp_Z.pkl","rb") as file:
            gan_amps = pickle.load(file)
            real_amps = pickle.load(file)
            z   = pickle.load(file)
            z2  = pickle.load(file)
            r   = pickle.load(file)
            r2  = pickle.load(file)
            phi = pickle.load(file)
            phi2= pickle.load(file)
            gan_se = pickle.load(file)
            real_se= pickle.load(file)
        file.close()
        print gan_amps.shape,real_amps.shape,r.shape,r2.shape
        print 'real amplitude:',[np.sum(real_amps[i]) for i in range(5)]
        print 'real raw energy:',real_se[:5]
        print 'GAN amplitude:',[np.sum(gan_amps[i]) for i in range(5)]
        print 'GAN raw energy:',gan_se[:5]
        #outliers = plot_amp_minus_energy(gan_se,real_se,gan_amps,real_amps)
        #plot_se_amp(gan_se,gan_amps,name='GAN')
        #plot_se_amp(real_se,real_amps,name='Real')
        #RE.plot([],gan_se,[],real_se,3500,9000)
        res_r,res_g,mu_r,mu_g = fit_gaussian(gan_se,real_se,[4000,7000,500], [6200,7500],[4000,8000,500],[6200+300,7500+300],True,FileName,unit='A.U.')
        #res_r,res_g,mu_r,mu_g = RE.fit_res(real_se,gan_se,[4000,7000,500,100], [4500,9000],[4000,7000,500,100],[4500,9000],True,FileName,unit='A.U.')
        #resfit_r,resfit_g,er,eg = RE.fit_res(real_amps,gan_amps,[4000,6500,500,3500], [5800,7400],[4000,6500,500,3500],[5800,7400],True,FileName,unit='A.U.')
        '''
        plt.rcParams['figure.figsize'] = [12, 8]
        fig, axs = plt.subplots(1,3)
        axs[0].hist(r[outliers],bins=40,label='R')
        axs[1].hist(phi[outliers],bins=40,label=r'$\phi$')
        axs[2].hist(z[outliers],bins=40,label='Z')
        axs[0].legend()
        axs[1].legend()
        axs[2].legend()
        plt.show()
        '''
        
        for i in range(0,3):
            plot_apd_z(gan_amps,real_amps,z,z2,i)
            plot_apd_r(gan_amps,real_amps,r,r2,i)
            plot_apd_phi(gan_amps,real_amps,phi,phi2,i)
        
    else:
        import ROOT
        ROOT.gSystem.Load("libEXOUtilities")
        ROOT.gSystem.Load("libEXOFitting")
        read(RunNumber,FileName,Phase)