##############################################################
#
# python ReconEvents.py --FileName flatten2 --ReadFile True/False --Phase 1/2
#
##############################################################

import os, re, sys, glob, math, random, pickle,csv
import numpy as np
import ROOT
import argparse
import matplotlib as mpl 
#mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats, special
from scipy import interpolate 
from scipy.optimize import curve_fit
from scipy.stats import kde
import random

def get_errs(n):
    err_list = []
    for i in range(0,len(n)):
        if n[i]!=0:
            err_i = np.sqrt(n[i])
        else:
            err_i = 1.
        err_list.append(err_i)
    return np.asarray(err_list)

def get_centers(bins,step):
    centers_list = []
    n = len(bins)-1
    for i in range(0,n):
        e_i = bins[i]+step/2.
        centers_list.append(e_i)
    return np.asarray(centers_list)

def read(runNumbers,filename,phase):
    # Under the directory: # gan_lightmap_recon/test_file_lightmap recon_with_lightmap/lightmap_2017_0nu-base64_v2
    recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/%s'%(phase,filename)
    #recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/%s/over2000'%(phase,filename)
    #recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/%s-lightmap-MPD-NoPooling-phase2-flatten-random'%(phase,filename)
    #recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/%s/before-moving-baseline/'%(phase,filename)
    #recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/gan_lightmap_recon/MPD-phase2-flatten-AmpPenalty3'%(phase)
    #recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase2/recon_with_lightmap/lightmap_2017_0nu-base64_v2/'
    #recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase2/gan_lightmap_recon/MPD-NoPooling-phase2-flatten/'
    t = ROOT.TChain('tree','tree')
    for runNumber in runNumbers:
        #t.Add(runNumber)
        print('Add file: %s'%runNumber)
        t.Add(recon_path+"/recon_%s*_gan_%s.root" % (runNumber,filename)) #recon_with_lightmap/
    ed = ROOT.EXOEventData()
    t.SetBranchAddress('EventBranch',ed)

    t2 = ROOT.TChain('tree','tree')
    for runNumber in runNumbers:
        #t2.Add(runNumber.replace('_gan_%s'%filename,'_real'))
        print('Add file: %s'%runNumber)
        t2.Add(recon_path+"/recon_%s*_real.root" % runNumber) #recon_with_lightmap/
    ed2 = ROOT.EXOEventData()
    t2.SetBranchAddress('EventBranch',ed2)
    print "Same nEvents?", t.GetEntries() == t2.GetEntries()

    gan_e_charge,gan_e_scint = [],[]
    wrong_event = []
    real_e_charge,real_e_scint = [],[]
    wrong_event_2 = []
    z,z2 = [],[]
    r,r2 = [],[]
    Id,Id2 = [],[]
    gan_amp,real_amp = [],[]

    for i in range(t.GetEntries()):
        t.GetEntry(i)
        t2.GetEntry(i)
        nev = ed.fEventNumber
        nev2 = ed2.fEventNumber

        ncc = ed.GetNumChargeClusters()
        ncc2 = ed2.GetNumChargeClusters()
        nsc = ed.GetNumScintillationClusters()
        nsc2 = ed2.GetNumScintillationClusters()
        if ncc != 1 or ncc2 != 1 or nsc != 1 or nsc2 != 1:
            wrong_event_2.append([i,nev,nev2,ncc,ncc2])
            wrong_event.append([i,nev,nev2,nsc,nsc2])
            continue
        ecc = 0
        for k in range(ncc):
            cc = ed.GetChargeCluster(k)
            #ecc += cc.fRawEnergy
            ecc += cc.fPurityCorrectedEnergy
            #ecc += cc.fCorrectedEnergy
            z.append(cc.fZ)
            r.append([cc.fX, cc.fY])
        gan_e_charge.append(ecc)
        ecc = 0
        for k in range(ncc2):
            cc = ed2.GetChargeCluster(k)
            #ecc += cc.fRawEnergy
            ecc += cc.fPurityCorrectedEnergy
            #ecc += cc.fCorrectedEnergy
            z2.append(cc.fZ)
            r2.append([cc.fX, cc.fY])
        real_e_charge.append(ecc)

        #if nsc != 1 or nsc2 != 1:
        #    wrong_event.append([i,nev,nev2,nsc,nsc2])
        #    continue
        esc = 0
        plane1,plane2 = 0,0
        for j in range(nsc):
            sc = ed.GetScintillationCluster(j)
            esc += sc.fRawEnergy
            wv1 = sc.GetAPDSignal(ROOT.EXOAPDSignal.kPlaneFit,1)
            wv2 = sc.GetAPDSignal(ROOT.EXOAPDSignal.kPlaneFit,2)
            
            if wv1: plane1 = wv1.fRawCounts
            if wv2: plane2 = wv2.fRawCounts
        gan_e_scint.append(esc)
        gan_amp.append(plane1+plane2)
        esc = 0
        plane1,plane2 = 0,0
        for j in range(nsc2):
            sc = ed2.GetScintillationCluster(j)
            esc += sc.fRawEnergy
            wv1 = sc.GetAPDSignal(ROOT.EXOAPDSignal.kPlaneFit,1)
            wv2 = sc.GetAPDSignal(ROOT.EXOAPDSignal.kPlaneFit,2)
            #plane1,plane2 = 0,0
            if wv1: plane1 = wv1.fRawCounts
            if wv2: plane2 = wv2.fRawCounts
        if esc == 6756.9145724690461 and filename == 'smallPhase1se':
            continue
        real_e_scint.append(esc)
        real_amp.append(plane1+plane2)

        Id.append(i)
        Id2.append(i)

    with open("./pkls/gan_real_E.pkl","wb") as file:
        pickle.dump(np.array(gan_e_charge),file)
        pickle.dump(np.array(gan_e_scint),file)
        pickle.dump(np.array(real_e_charge),file)
        pickle.dump(np.array(real_e_scint),file)
        pickle.dump(np.array(wrong_event),file)
        pickle.dump(np.array(wrong_event_2),file)
        pickle.dump(np.array(z),file)
        pickle.dump(np.array(z2),file)
        pickle.dump(np.array(Id),file)
        pickle.dump(np.array(Id2),file)
        pickle.dump(np.array(r),file)
        pickle.dump(np.array(r2),file)
        pickle.dump(np.array(gan_amp),file)
        pickle.dump(np.array(real_amp),file)
    
    print "wrong_event",wrong_event[:20]
    print "wrong_event_r",wrong_event_2[:20]
    print "gan_e_charge",len(gan_e_charge)
    print "gan_e_scint",len(gan_e_scint)
    print "real_e_charge",len(real_e_charge)
    print "real_e_scint",len(real_e_scint)
    
    return gan_e_charge,gan_e_scint,real_e_charge,real_e_scint
def plot(gan_e_charge,gan_e_scint,real_e_charge,real_e_scint,charge_range,scint_range):
    left = 100
    right = charge_range
    nbins = (right - left)/ 40
    mybins = np.linspace(left,right,nbins+1)
    e_step = (right-left)/nbins
    nGANeCharge, binsGANeCharge = np.histogram(gan_e_charge,bins=mybins)
    nRealeCharge, binsRealeCharge = np.histogram(real_e_charge,bins=mybins)

    err_gec = get_errs(nGANeCharge)
    err_rec = get_errs(nRealeCharge)
    cent  = get_centers(mybins,e_step)
    plt.rcParams['figure.figsize'] = [12, 10]
    plt.errorbar(cent,nGANeCharge,yerr=err_gec,ecolor='lightgrey',fmt='o',color='r',markersize=6.,label='GAN charge E')
    plt.errorbar(cent,nRealeCharge,yerr=err_rec,ecolor='lightgrey',fmt='o',color='b',markersize=6.,label='Real charge E')
    plt.legend()
    plt.title(r'Charge E',fontsize=20)
    plt.xlabel(r'fCorrectedEnergy [keV]',fontsize=20)
    plt.ylabel(r'Counts/(%.0f keV)'%(e_step),fontsize=20)
    plt.savefig('pics/charge_E_gan_real.png',bbox_inches = 'tight')
    plt.show()

    left = 100
    right = scint_range
    nbins = (right - left) / 40
    mybins = np.linspace(left,right,nbins+1)
    e_step = (right-left)/nbins
    nGANeScint, binsGANeScint = np.histogram(gan_e_scint,bins=mybins)
    nRealeScint, binsRealeScint = np.histogram (real_e_scint,bins=mybins)

    err_ges = get_errs(nGANeScint)
    err_res = get_errs(nRealeScint)
    cent  = get_centers(mybins,e_step)
    plt.rcParams['figure.figsize'] = [12, 10]
    plt.errorbar(cent,nGANeScint,yerr=err_ges,ecolor='lightgrey',fmt='o',color='r',markersize=6.,label='GAN light E')
    plt.errorbar(cent,nRealeScint,yerr=err_res,ecolor='lightgrey',fmt='o',color='b',markersize=6.,label='Real light E')
    plt.legend()
    plt.title(r'light E',fontsize=20)
    plt.xlabel(r'fRawEnergy [A.U.]',fontsize=20)
    plt.ylabel(r'Counts/(%.0f A.U.)'%(e_step),fontsize=20)
    plt.savefig('pics/light_E_gan_real.png',bbox_inches = 'tight')
    plt.show()


    plt.scatter(gan_e_charge,real_e_charge,color='r',label = "charge")
    plt.legend()
    plt.title("fPurityCorrectedEnergy")
    plt.xlabel(r"GAN E [keV]",fontsize=20)
    plt.ylabel(r"Real E [keV]",fontsize=20)
    plt.savefig("pics/event_chargeE_gan_real.png",bbox_inches = 'tight')
    plt.show()


    plt.scatter(gan_e_scint,real_e_scint,color='b',label = "light")
    plt.xlim([0,10000])
    plt.ylim([0,10000])
    plt.legend()
    plt.xlabel(r"GAN E [keV]",fontsize=20)
    plt.ylabel(r"Real E [keV]",fontsize=20)
    plt.savefig("pics/event_lightE_gan_real.png",bbox_inches = 'tight')
    plt.show()
def plot_corr(gan_e_charge,gan_e_scint,real_e_charge,real_e_scint):
    plt.rcParams['figure.figsize'] = [12, 10]
    nr,ng = real_e_scint.shape[0],gan_e_scint.shape[0]
    xmax,ymax = gan_e_charge.max(),gan_e_scint.max()#1500,2000 #
    xmin,ymin = 2000,2000#gan_e_charge.min(),gan_e_scint.min()
    
    plt.scatter(real_e_charge[:10000],real_e_scint[:10000],color='b',label = "real")
    plt.scatter(gan_e_charge[:10000],gan_e_scint[:10000],color='r',label = "GAN")
    plt.legend()
    #plt.title("fPurityCorrectedEnergy")
    plt.xlabel(r"Charge E [keV]",fontsize=20)
    plt.ylabel(r"Light E [keV]",fontsize=20)
    plt.savefig("pics/gan_corr.png",bbox_inches = 'tight')
    plt.show()
    

    fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(11, 5))
    nbins = 50
    print(real_e_charge.shape,real_e_scint.shape)
    data = np.concatenate((real_e_charge.reshape(nr,1),real_e_scint.reshape(nr,1)),axis=1)
    np.random.shuffle(data)
    x,y = data.T[:40000]
    k = kde.gaussian_kde([x,y])
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    axes[0].set_title('Real')
    axes[0].pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.BuGn_r)
    axes[0].contour(xi, yi, zi.reshape(xi.shape) )
    axes[0].set_xlabel('Charge E [keV]')
    axes[0].set_ylabel('Light E [keV]')

    data = np.concatenate((gan_e_charge.reshape(ng,1),gan_e_scint.reshape(ng,1)),axis=1)
    np.random.shuffle(data)
    x,y = data.T[:40000]
    k = kde.gaussian_kde([x,y])
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    axes[1].set_title('GAN')
    axes[1].pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gouraud', cmap=plt.cm.BuGn_r)
    axes[1].contour(xi, yi, zi.reshape(xi.shape) )
    axes[1].set_xlabel('Charge E [keV]')
    plt.savefig("pics/gan_contour.png",bbox_inches = 'tight')
    plt.show()

def fit_2D_Co_peaks(rce,rse,gce,gse,FileName):
    nr,ng = rse.shape[0],gse.shape[0]
    xmax,ymax = 1400,6000 #gan_e_charge.max(),gan_e_scint.max()
    xmin,ymin = 900,3000
    nbins = 201
    # Fit for GAN
    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(11, 5))
    # create data
    data = np.concatenate((gce.reshape(nr,1),gse.reshape(nr,1)),axis=1)
    np.random.shuffle(data)
    data = np.array([i for i in data if i[1] > -8*(i[0]-1000) +3200 and i[1] < -8*(i[0]-1197.58) + 3906.25])
    x,y = data.T[:2000]
    k = kde.gaussian_kde([x,y])
    # Create x and y indices
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    # add some noise to the data and try to fit the data generated beforehand
    initial_guess = [  3.92519766e-06,   1.10162384e+03,   4.01306343e+03,   5.23940438e+02,
   7.77859520e+01,   1.52376229e+00,  -6.58313785e-09]
    popt, pcov = curve_fit(twoD_Gaussian, (xi, yi), zi, p0=initial_guess)
    print(popt)
    gan_peak1 = popt[2]

    data_fitted = twoD_Gaussian((xi, yi), *popt)
    #fig, ax = plt.subplots(1, 1)
    ax[0].hold(True)
    ax[0].pcolormesh(xi, yi, zi.reshape(201,201))
    ax[0].scatter([popt[1]],[popt[2]],marker = '*')
    ax[0].annotate("(%.2f,%.2f)"%(popt[1],popt[2]), (popt[1],popt[2]),color=(202. / 255., 121. / 256., 0. / 255.))
    #ax.imshow(zi.reshape(201, 201), cmap=plt.cm.jet, origin='bottom',
    #    extent=(xmin, xmax, ymin, ymax))
    ax[0].contour(xi, yi, data_fitted.reshape(201, 201), 8, colors='w')

    # create data
    data = np.concatenate((gce.reshape(nr,1),gse.reshape(nr,1)),axis=1)
    np.random.shuffle(data)
    data = np.array([i for i in data if i[1] > -8*(i[0]-1197.58) + 3906.25])
    x,y = data.T[:2000]
    k = kde.gaussian_kde([x,y])
    # Create x and y indices
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    # add some noise to the data and try to fit the data generated beforehand
    initial_guess = [  4.57803849e-06,   1.27062772e+03,   4.10203847e+03,   5.68545770e+02,
   6.05877185e+01,   1.49204124e+00,   1.32866751e-08]
    popt, pcov = curve_fit(twoD_Gaussian, (xi, yi), zi, p0=initial_guess)
    print(popt)
    gan_peak2 = popt[2]

    data_fitted = twoD_Gaussian((xi, yi), *popt)
    #fig, ax = plt.subplots(1, 1)
    ax[1].hold(True)
    ax[1].pcolormesh(xi, yi, zi.reshape(201,201))
    ax[1].scatter([popt[1]],[popt[2]],marker = '*')
    ax[1].annotate("(%.2f,%.2f)"%(popt[1],popt[2]), (popt[1],popt[2]),color=(202. / 255., 121. / 256., 0. / 255.))
    #ax.imshow(zi.reshape(201, 201), cmap=plt.cm.jet, origin='bottom',
    #    extent=(xmin, xmax, ymin, ymax))
    ax[1].contour(xi, yi, data_fitted.reshape(201, 201), 8, colors='w')
    fig.suptitle('GAN Co Peaks')
    plt.savefig('pics/Co_contour_fit_gan.png')
    plt.show()
    #####################################################################
    # Fit for real 
    fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(11, 5))
    # create data
    data = np.concatenate((rce.reshape(nr,1),rse.reshape(nr,1)),axis=1)
    np.random.shuffle(data)
    data = np.array([i for i in data if i[1] > -8*(i[0]-1000) +3200 and i[1] < -8*(i[0]-1198.39) + 4350])
    x,y = data.T[:2000]
    k = kde.gaussian_kde([x,y])
    # Create x and y indices
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    # add some noise to the data and try to fit the data generated beforehand
    initial_guess = [  3.77390767e-06,   1.10343077e+03,   4.02100232e+03,   5.29456482e+02,
   8.05883776e+01,   1.52869298e+00,  -9.35707928e-09]
    popt, pcov = curve_fit(twoD_Gaussian, (xi, yi), zi, p0=initial_guess)
    print(popt)
    real_peak1 = popt[2]

    data_fitted = twoD_Gaussian((xi, yi), *popt)
    #fig, ax = plt.subplots(1, 1)
    ax[0].hold(True)
    ax[0].pcolormesh(xi, yi, zi.reshape(201,201))
    ax[0].scatter([popt[1]],[popt[2]],marker = '*')
    ax[0].annotate("(%.2f,%.2f)"%(popt[1],popt[2]), (popt[1],popt[2]),color=(202. / 255., 121. / 256., 0. / 255.))
    #ax.imshow(zi.reshape(201, 201), cmap=plt.cm.jet, origin='bottom',
    #    extent=(xmin, xmax, ymin, ymax))
    ax[0].contour(xi, yi, data_fitted.reshape(201, 201), 8, colors='w')

    # create data
    data = np.concatenate((rce.reshape(nr,1),rse.reshape(nr,1)),axis=1)
    np.random.shuffle(data)
    data = np.array([i for i in data if i[1] > -8*(i[0]-1198.39) + 4350])
    x,y = data.T[:2000]
    k = kde.gaussian_kde([x,y])
    # Create x and y indices
    xi, yi = np.mgrid[xmin:xmax:nbins*1j, ymin:ymax:nbins*1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))

    # add some noise to the data and try to fit the data generated beforehand
    initial_guess = [  5.17231495e-06,   1.27335375e+03,   4.55431617e+03,   5.29730653e+02,
   5.93918563e+01,   1.49070174e+00,  -1.25562248e-08]
    popt, pcov = curve_fit(twoD_Gaussian, (xi, yi), zi, p0=initial_guess)
    print(popt)
    real_peak2 = popt[2]

    data_fitted = twoD_Gaussian((xi, yi), *popt)
    #fig, ax = plt.subplots(1, 1)
    ax[1].hold(True)
    ax[1].pcolormesh(xi, yi, zi.reshape(201,201))
    ax[1].scatter([popt[1]],[popt[2]],marker = '*')
    ax[1].annotate("(%.2f,%.2f)"%(popt[1],popt[2]), (popt[1],popt[2]),color=(202. / 255., 121. / 256., 0. / 255.))
    #ax.imshow(zi.reshape(201, 201), cmap=plt.cm.jet, origin='bottom',
    #    extent=(xmin, xmax, ymin, ymax))
    ax[1].contour(xi, yi, data_fitted.reshape(201, 201), 8, colors='w')
    fig.suptitle('Real Co Peaks')
    plt.savefig('pics/Co_contour_fit_real.png')
    plt.show()

    '''
    # plot twoD_Gaussian data generated above
    plt.figure()
    #
    X= np.linspace(600,xmax,201)
    Y = -8*(X-1000) +3200
    Y2 = -8*(X-1205.65) + 4343.75
    Y3 = -8*(X-1197.58) + 3906.25
    plt.plot(X,Y,color='r')
    plt.plot(X,Y2,color='b')
    plt.pcolormesh(xi, yi, zi.reshape(201,201))
    plt.colorbar()
    plt.savefig('pics/Co_contour.png')
    plt.clf()
    '''
    return gan_peak1,gan_peak2,real_peak1,real_peak2


def cut(arr,arr2):
    selected_arr,selected_arr2 = [],[]
    index_list = []
    n = min(len(arr),len(arr2))
    for i in range(n):
        if arr[i] > 1.1 * arr2[i]:
            selected_arr.append(arr[i])
            selected_arr2.append(arr2[i])
            index_list.append(i)
    with open("./pkls/gan_real_cut_index.pkl","wb") as file:
        pickle.dump(index_list,file)
    file.close()
    return selected_arr,selected_arr2,index_list

def plot_z(z,z2):
    left = -210
    right = 210
    nbins = 40
    mybins = np.linspace(left,right,nbins+1)
    e_step = (right-left)/nbins
    nGANeCharge, binsGANeCharge = np.histogram(z,bins=mybins)
    nRealeCharge, binsRealeCharge = np.histogram(z2,bins=mybins)

    err_gec = get_errs(nGANeCharge)
    err_rec = get_errs(nRealeCharge)
    cent  = get_centers(mybins,e_step)
    plt.rcParams['figure.figsize'] = [12, 10]
    plt.errorbar(cent,nGANeCharge,yerr=err_gec,ecolor='lightgrey',fmt='o',color='r',markersize=6.,label='GAN charge Z')
    plt.errorbar(cent,nRealeCharge,yerr=err_rec,ecolor='lightgrey',fmt='o',color='b',markersize=6.,label='Real charge Z')
    plt.legend()
    plt.title(r'Charge Z',fontsize=20)
    plt.xlabel(r'fZ [cm]',fontsize=20)
    plt.ylabel(r'Counts/(%.0f cm)'%(e_step),fontsize=20)
    plt.savefig('pics/charge_Z_gan_real_cut.png',bbox_inches = 'tight')
    plt.show()
#get chi-square and number of degrees of freedom
def my_chisq(func,popt,x,data,data_err,cutoff_l=0,cutoff_h=np.inf,pearson=False):
    chi=0
    ndf = -1.*len(popt)
    
    for i in range(0,len(data)):
        if x[i]<cutoff_l or x[i]>cutoff_h:
            continue
        diff = (func(x,*popt)[i]-data[i])**2
        err = func(x,*popt)[i] if pearson==True else data_err[i]**2
        chi = chi + diff/err
        #print diff, err, chi
        ndf = ndf + 1
    return chi, ndf
#make error bars and bin-centered points
def get_errs(n):
    err_list = []
    for i in range(0,len(n)):
        if n[i]!=0:
            err_i = np.sqrt(n[i])
        else:
            err_i = 1.
        err_list.append(err_i)
    return np.asarray(err_list)

def get_centers(bins,step):
    centers_list = []
    n = len(bins)-1
    for i in range(0,n):
        e_i = bins[i]+step/2.
        centers_list.append(e_i)
    return np.asarray(centers_list)
#fit function
'''
def g_lin_fit(x,a,m,s,p0):
        #super inefficent way to use python:
        y = []
        for x_i in x:
            y_i= a/s/np.sqrt(2*np.pi)*np.exp(-(x_i-m)*(x_i-m)/2/s/s) + math.erf((m - x_i)/s)*p0 + p0
            y.append(y_i)
        return y
'''
def g_lin_fit(x,a,m,s,p0):
    return a/s/np.sqrt(2*np.pi)*np.exp(-(x-m)*(x-m)/2/s/s) + special.erfc((x - 0.91*m)/s)*p0 
def gaussian(x,a,m,s):
    return a/s/np.sqrt(2*np.pi)*np.exp(-(x-m)*(x-m)/2/s/s)
def Erfc(x,m,s,p0):
    return special.erfc((x - 0.91*m)/s)*p0
def double_peak(x,a1,m1,s1,a2,m2,s2):
    return a1/s1/np.sqrt(2*np.pi)*np.exp(-(x-m1)*(x-m1)/2/s1/s1) + a2/s2/np.sqrt(2*np.pi)*np.exp(-(x-m2)*(x-m2)/2/s2/s2)
#define model function and pass independant variables x and y as a list
def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

def double_two_Gaussian((x,y),amplitude, xo, yo, sigma_x, sigma_y, theta, offset,amplitude1, x1, y1, sigma_x1, sigma_y1, theta1, offset1):
    return twoD_Gaussian((x,y),amplitude, xo, yo, sigma_x, sigma_y, theta, offset) + twoD_Gaussian((x,y),amplitude1, x1, y1, sigma_x1, sigma_y1, theta1, offset1)
def poly2(x,a,b,c):
    y = []
    for xi in x:
        yi = a * (xi - b) ** 2 + c
        y.append(yi)
    return y
#get subset of array within the fit range
def get_range(cents,n,err,left,right):
    centv         = np.vstack((cents,n)).T
    cent_err      = np.vstack((cents,err)).T
    cent_range    = centv[((centv[:,0]>left)&(centv[:,0]<right)),0]
    n_range       = centv[((centv[:,0]>left)&(centv[:,0]<right)),1]
    err_range     = cent_err[((cent_err[:,0]>left)&(cent_err[:,0]<right)),1]
    return cent_range, n_range, err_range

def fit_res(reals,gens,first,bounds,gan_first,gan_bounds,show=False,name = '',unit='keV'):
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
    ngens, bins   = np.histogram(gens,  bins=mybins)
    
    err_r = get_errs(nreals)
    err_g = get_errs(ngens)
    cent  = get_centers(mybins,e_step)

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
    cent_range_r, nrange_r, err_range_r = get_range(cent,nreals,err_r,left_r,right_r)
    cent_range_g, nrange_g, err_range_g = get_range(cent,ngens,err_g,left_g,right_g)
    xr = np.linspace(left_r,right_r,500)
    xg = np.linspace(left_g,right_g,500)

    #initial parameter guess and bounds
    #first = (4000,8000,50,100)
    bound=([0,2400,20,0],[np.inf,2700,300,5000])

    #now fit and plot results
    popt_r, pcov = curve_fit(g_lin_fit, cent_range_r, nrange_r,first, maxfev=10000,sigma=err_range_r)
    perr_r = np.sqrt(np.diag(pcov))
    mychi2,ndof = my_chisq(g_lin_fit,popt_r,cent_range_r,nrange_r,err_range_r)
    plt.plot(xr, g_lin_fit(xr, *popt_r), 'b-',linewidth=5.0,label=r'$\mu$=%5.1f$\pm$%3.1f %s, $\sigma$/$\mu$=%1.1f$\pm$%1.1f%%,  $\chi$$^2$/ndf=%.1f' % (popt_r[1],perr_r[1],unit,popt_r[2]/popt_r[1]*100.,perr_r[2]/popt_r[1]*100.,mychi2/ndof) )
    plt.plot(xr, gaussian(xr, *popt_r[:3]), 'b--',linewidth=1.0)
    plt.plot(xr, Erfc(xr, *popt_r[1:]), 'b--',linewidth=1.0) 
    res_r = popt_r[2]/popt_r[1]*100.
    mu_r = popt_r[1]

    popt_g, pcov = curve_fit(g_lin_fit, cent_range_g, nrange_g,gan_first, maxfev=10000,sigma=err_range_g)
    perr_g = np.sqrt(np.diag(pcov))
    mychi2,ndof = my_chisq(g_lin_fit,popt_g,cent_range_g,nrange_g,err_range_g)
    plt.plot(xg, g_lin_fit(xg, *popt_g), 'r-',linewidth=5.0,label=r'$\mu$=%5.1f$\pm$%3.1f %s, $\sigma$/$\mu$=%1.1f$\pm$%1.1f%%,  $\chi$$^2$/ndf=%.1f' % (popt_g[1],perr_g[1],unit,popt_g[2]/popt_g[1]*100.,perr_g[2]/popt_g[1]*100.,mychi2/ndof) )
    plt.plot(xr, gaussian(xr, *popt_g[:3]), 'r--',linewidth=1.0 )
    plt.plot(xr, Erfc(xr, *popt_g[1:]), 'r--',linewidth=1.0 )
    res_g = popt_g[2]/popt_g[1]*100.
    mu_g = popt_g[1]

    plt.legend(loc='upper left',numpoints = 1,frameon=True,prop={'size': 16})
    plt.gca().tick_params(labelsize=18)
    #plt.title('light E spectrum for %s'%name)

    plt.savefig('./pics/reco_fit_check_%s.png'%name,bbox_inches = 'tight')
    if show:
        plt.show()
    else:
        plt.clf()
    return res_r,res_g,mu_r,mu_g
def rot_res(rce,rse,gce,gse):
    theta = np.linspace(1. ,1.4,10)
    #initial parameter guess and bounds
    res_r,res_g = [],[]
    for i in theta:
        re = rce* np.sin(i) + rse * np.cos(i)
        ge = gce* np.sin(i) + gse * np.cos(i)
        bounds_real = [2200*np.cos(i)+2200*np.sin(i),2800*np.cos(i)+2800*np.sin(i)]
        first_real = [4000,2615*np.sin(i)+2615*np.cos(i),50,3500] # p3 + erf((p1 - x)/p2) * p3 + p0/p1/sqrt(2Pi) * exp(-(x-p1)^2 / (2*p2^2))
        bounds_gan = [2200*np.cos(i)+2200*np.sin(i),2800*np.cos(i)+2800*np.sin(i)]
        first_gan = [4000,2615*np.sin(i)+2615*np.cos(i),50,3500] # p3 + erf((p1 - x)/p2) * p3 + p0/p1/sqrt(2Pi) * exp(-(x-p1)^2 / (2*p2^2))
        shift = [0,0]
        resfit_r,resfit_g,_,_ = fit_res(re,ge,first_real,bounds_real,first_gan,bounds_gan,False)
        res_r.append(resfit_r)
        res_g.append(resfit_g)

    plt.scatter(theta,res_r,label = "Real Resolution",color = 'b')
    plt.scatter(theta,res_g,label = "GAN Resolution",color = 'r')

    first = [5,1.3,2.5]
    popt_r, pcov = curve_fit(poly2, theta, np.array(res_r),first, maxfev=10000)
    popt_g, pcov = curve_fit(poly2, theta, np.array(res_g),first, maxfev=10000)
    plt.plot(theta, poly2(theta, *popt_r), 'b-',label=r'$\theta_{min}$=%.3f' %(popt_r[1]))
    plt.plot(theta, poly2(theta, *popt_g), 'r-',label=r'$\theta_{min}$=%.3f' %(popt_g[1]))

    plt.ylabel(r"$\sigma$/$\mu$ $\%$")
    plt.xlabel(r"rotation angle $\theta$")
    plt.legend()
    plt.savefig("./pics/rotation_fit_check.png",bbox_inches = 'tight')
    plt.show()
    return popt_r[1],popt_g[1]
def cali_res(rse,gse,real_first,real_bounds,gan_first,gan_bounds,name):
    # bounds for real histogram; bounds - shift for gan histogram
    #bounds = [6000,10000]
    #first = [4000,8000,50,100]
    #shift = [1000,1000]
    resfit_r,resfit_g,er,eg = fit_res(rse,gse,real_first,real_bounds,gan_first,gan_bounds,True,name)

    with open("./pkls/gan_real_cali_E.pkl","wb") as file:
        pickle.dump(gse / eg * 2615,file)
        pickle.dump(rse / er * 2615,file)
    return resfit_r,resfit_g,er,eg
#
def fit_Co_peak(reals,gens,first,bounds,gan_first,gan_bounds,show=False,name = ''):
    #set the figure size
    plt.rcParams['figure.figsize'] = [12, 10]

    #create binning
    left = 0.
    right = bounds[1] + 200#12000
    nbins = 200#(right - left) / 20
    mybins = np.linspace(left,right,nbins+1)
    e_step = (right-left)/nbins

    #get bin content
    nreals, bins  = np.histogram(reals, bins=mybins)
    ngens, bins   = np.histogram(gens,  bins=mybins)
    
    err_r = get_errs(nreals)
    err_g = get_errs(ngens)
    cent  = get_centers(mybins,e_step)

    #now actually plot
    #plt.step(cent,nlabels,'k',where='mid',linewidth=2) #try adding this for coarse binning - may look nicer
    plt.errorbar(cent,nreals,yerr=err_r,ecolor='lightgrey',capsize=0,elinewidth=2,fmt='o',color='b',markersize=4.,label=r'Light E (Reco on Real Waveforms)')
    plt.errorbar(cent,ngens,yerr=err_g,ecolor='lightgrey',capsize=0,elinewidth=2,fmt='o',color='r',markersize=4.,label=r'Light E (Reco on GAN Waveforms)')

    #Accuracy of axis labels tells people if you are a professional or a hack!
    plt.xlabel(r'E (keV)',fontsize=20)
    plt.ylabel(r'Counts/(%.0f keV)'%(e_step),fontsize=20)

    #Set axes limits
    plt.axis((left,right,0,np.max(nreals)*1.75))

    #fit range
    left_r  = bounds[0]
    right_r = bounds[1]
    left_g  = gan_bounds[0]
    right_g = gan_bounds[1]
    cent_range_r, nrange_r, err_range_r = get_range(cent,nreals,err_r,left_r,right_r)
    cent_range_g, nrange_g, err_range_g = get_range(cent,ngens,err_g,left_g,right_g)
    xr = np.linspace(left_r,right_r,500)
    xg = np.linspace(left_g,right_g,500)

    #now fit and plot results
    popt_r, pcov = curve_fit(double_peak, cent_range_r, nrange_r,first, maxfev=10000,sigma=err_range_r)
    perr_r = np.sqrt(np.diag(pcov))
    mychi2,ndof = my_chisq(double_peak,popt_r,cent_range_r,nrange_r,err_range_r)
    plt.plot(xr, double_peak(xr, *popt_r), 'b-',label=r'''$\mu_1$=%5.1f$\pm_1$%3.1f keV, $\sigma_1$/$\mu_1$=%1.1f$\pm$%1.1f%%,$\chi$$^2$/ndf=%.1f 
    $\mu_2$=%5.1f$\pm_2$%3.1f keV, $\sigma_2$/$\mu_2$=%1.1f$\pm$%1.1f%%,  $\chi$$^2$/ndf=%.1f  ''' 
                                                           %(popt_r[1],perr_r[1],popt_r[2]/popt_r[1]*100.,perr_r[2]/popt_r[1]*100.,mychi2/ndof,
                                                             popt_r[4],perr_r[4],popt_r[5]/popt_r[4]*100.,perr_r[5]/popt_r[1]*100.,mychi2/ndof) )
    res_r1,res_r2 = popt_r[2]/popt_r[1]*100.,popt_r[5]/popt_r[4]*100.
    mu_r1,mu_r2 = popt_r[1],popt_r[4]

    popt_g, pcov = curve_fit(double_peak, cent_range_g, nrange_g,gan_first, maxfev=10000,sigma=err_range_g)
    perr_g = np.sqrt(np.diag(pcov))
    mychi2,ndof = my_chisq(double_peak,popt_g,cent_range_g,nrange_g,err_range_g)
    plt.plot(xg, double_peak(xg, *popt_g), 'r-',label=r'''$\mu_1$=%5.1f$\pm_1$%3.1f keV, $\sigma_1$/$\mu_1$=%1.1f$\pm$%1.1f%%,$\chi$$^2$/ndf=%.1f 
    $\mu_2$=%5.1f$\pm_2$%3.1f keV, $\sigma_2$/$\mu_2$=%1.1f$\pm$%1.1f%%,  $\chi$$^2$/ndf=%.1f  ''' 
                                                           %(popt_g[1],perr_g[1],popt_g[2]/popt_g[1]*100.,perr_g[2]/popt_g[1]*100.,mychi2/ndof,
                                                             popt_g[4],perr_g[4],popt_g[5]/popt_g[4]*100.,perr_g[5]/popt_g[1]*100.,mychi2/ndof) )
    res_g1,res_g2 = popt_g[2]/popt_g[1]*100.,popt_g[5]/popt_g[4]*100.,
    mu_g1,mu_g2 = popt_g[1],popt_g[4]

    plt.legend(loc='upper left',numpoints = 1,frameon=True,prop={'size': 16})
    plt.gca().tick_params(labelsize=18)
    plt.title('light E spectrum for %s'%name)

    plt.savefig('./pics/reco_fit_check_%s.png'%name,bbox_inches = 'tight')
    if show:
        plt.show()
    else:
        plt.clf()
    return res_r1,res_r2,res_g1,res_g2,mu_r1,mu_r2,mu_g1,mu_g2

def rot_Co_res(rce,rse,gce,gse,FileName):
    theta = np.linspace(np.pi * 0.22 ,np.pi * 0.27,11)
    #initial parameter guess and bounds
    res_r1,res_r2,res_g1,res_g2 = [],[],[],[]
    for i in theta:
        re = rce* np.sin(i) + rse * np.cos(i)
        ge = gce* np.sin(i) + gse * np.cos(i)
        #bounds = [2000*np.cos(i)+2000*np.sin(i),3500*np.cos(i)+3500*np.sin(i)]
        #first = (4000,2615*np.cos(i) + 2615 * np.sin(i),50,100) # p3 + erf((p1 - x)/p2) * p3 + p0/p1/sqrt(2Pi) * exp(-(x-p1)^2 / (2*p2^2))
        first_real = [4000,4000*np.sin(i) + 1180*np.cos(i),100,4000,4500*np.sin(i)+1300*np.cos(i),100]
        bounds_real = [1000*np.sin(i)+500*np.cos(i),8000*np.sin(i)+1500*np.cos(i)]
        first_gan = [4000,3000*np.sin(i)+1180*np.cos(i),100,4000,4500*np.sin(i)+1300*np.cos(i),100]
        bounds_gan = [1000*np.sin(i)+500*np.cos(i),8000*np.sin(i)+1500*np.cos(i)]
        resfit_r1,resfit_r2,resfit_g1,resfit_g2,peak_real1,peak_real2,peak_gan1,peak_gan1 = fit_Co_peak(re,ge,first_real,bounds_real,first_gan,bounds_gan,True,FileName)
        print(r'Angle %.2f $\pi$'%(i/np.pi),peak_real1,peak_real2,peak_gan1,peak_gan1)
        res_r1.append(resfit_r1)
        res_r2.append(resfit_r2)
        res_g1.append(resfit_g1)
        res_g2.append(resfit_g2)

    plt.scatter(theta,res_r1,label = "Real Resolution 1",color = 'b',marker = 'o')
    plt.scatter(theta,res_r2,label = "Real Resolution 2",color = 'b',marker = 'x')
    plt.scatter(theta,res_g1,label = "GAN Resolution 1",color = 'r',marker = 'o')
    plt.scatter(theta,res_g2,label = "GAN Resolution 2",color = 'r',marker = 'x')

    first = [5,1.1,1.5]
    '''
    popt_r1, pcov = curve_fit(poly2, theta, np.array(res_r1),first, maxfev=10000)
    popt_r2, pcov = curve_fit(poly2, theta, np.array(res_r2),first, maxfev=10000)
    popt_g1, pcov = curve_fit(poly2, theta, np.array(res_g1),first, maxfev=10000)
    popt_g2, pcov = curve_fit(poly2, theta, np.array(res_g2),first, maxfev=10000)
    plt.plot(theta, poly2(theta, *popt_r1), 'b-',label=r'$\theta_{min}$=%.3f' %(popt_r1[1]))
    plt.plot(theta, poly2(theta, *popt_r2), 'b--',label=r'$\theta_{min}$=%.3f' %(popt_r2[1]))
    plt.plot(theta, poly2(theta, *popt_g1), 'r-',label=r'$\theta_{min}$=%.3f' %(popt_g1[1]))
    plt.plot(theta, poly2(theta, *popt_r2), 'r--',label=r'$\theta_{min}$=%.3f' %(popt_r2[1]))
    '''
    plt.ylabel(r"$\sigma$/$\mu$ $\%$")
    plt.xlabel(r"rotation angle $\theta$")
    plt.legend()
    plt.savefig("./pics/rotation_fit_check.png",bbox_inches = 'tight')
    plt.show()

    # Co-60 122keV and 136 keV

    with open("./pkls/gan_real_cali_rotated_E.pkl","wb") as file:
        pickle.dump(gse / eg * 2615,file)
        pickle.dump(rse / er * 2615,file)
    return resfit_r,resfit_g,er,eg

    return popt_r[1],popt_g[1]
def quadratic_function(E,a,b,c):    
    B = (a*(E**2.0)) + (b*E) + c
    return B
def quadratic_fit(Co_1,Co_2,Th):
    x = np.array([Co_1,Co_2,Th])
    y = np.array([1173,1332,2615])
    popt = np.polyfit(x, y, 2)
    print(popt)
    f = np.poly1d(popt)
    xi = np.linspace(500,9000,201)
    yi = f(xi)
    plt.plot(xi,yi,label=r'$E_c$ = %.2f $E_{raw}^2$ + %.2f $E_{raw}$ +%.2f'%(popt[0],popt[1],popt[2]))
    plt.scatter(x,y,marker = 'x')
    plt.legend()
    plt.show()
    return popt
def calibration_with_quadratic_fit(rse,gse,FileName):
    popt_gan = quadratic_fit(3.55005150e+03,4.10203813e+03,8022.1)
    popt_real = quadratic_fit(4.01747007e+03,4.55010048e+03,8618.5)
    f_gan = np.poly1d(popt_gan)
    f_real = np.poly1d(popt_real)
    with open("./pkls/gan_real_cali_E.pkl","wb") as file:
        pickle.dump(f_gan(gse),file)
        pickle.dump(f_real(rse),file)
def fiducial_cut(res,ges,r,r2,z,z2):
    n,n2 = len(res),len(ges)
    real_cut,gan_cut = [],[]
    # (r2[i][0]>120) and (r2[i][0]<150)  and np.abs(r2[i][1])<50 and (np.abs(z2[i]) < 50) and (np.abs(z2[i]) > 30)
    real_cut = [res[i] for i in range(n) if  (r2[i][0]**2 + r2[i][1]**2<168**2) and (np.abs(z2[i])<192.2)]
    gan_cut  = [ges[i] for i in range(n2) if (r2[i][0]**2 + r2[i][1]**2<168**2) and (np.abs(z2[i])<192.2)]
    return np.array(real_cut),np.array(gan_cut)
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #parser.add_argument('--RunNumber', type=int, nargs='+', default=5872)
    parser.add_argument('--FileName', type=str, nargs=1, default='unflatten')
    parser.add_argument('--ReadFile', type=bool,default=False)
    parser.add_argument('--Phase', type=int,default=1)
    args = parser.parse_args()
    #RunNumber = args.RunNumber
    FileName = args.FileName[0]
    readfile = args.ReadFile
    Phase = args.Phase
    print readfile

    RunNumber = []
    with open('phase%i_runset.csv'%Phase,'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if (row[1] == 'wTh') and row[3] == 'test' and row[0] > 3717:
                #print row[2]
                RunNumber.append(row[0])

    if readfile:
        with open("./pkls/gan_real_E.pkl","rb") as file:
            gec = pickle.load(file)
            ges = pickle.load(file)
            rec = pickle.load(file)
            res = pickle.load(file)
            w_s = pickle.load(file)
            w_c = pickle.load(file)
            z   = pickle.load(file)
            z2  = pickle.load(file)
            Id   = pickle.load(file)
            Id2  = pickle.load(file)
            r   = pickle.load(file)
            r2  = pickle.load(file)
            gan_amp = pickle.load(file)
            real_amp = pickle.load(file)
        file.close()
        print len(r),len(r2),len(ges),len(res)

        #rot_Co_res(rec,res,gec,ges,FileName)
        #gan_peak1,gan_peak2,real_peak1,real_peak2=fit_2D_Co_peaks(rec,res,gec,ges,FileName)
        #plot_corr(gec,ges,rec,res)
        #real_res1,real_res2,gan_res1,gan_res2,er1,er2,eg1,eg2 = fit_Co_peak(res,ges,[4000,3000,100,4000,4000,100],[1000,6000],[4000,3000,100,4000,3500,100],[1000,6000],FileName)

        ## For phase 1
        plot(gec,ges,rec,res,3500,9000)
        #res,ges = fiducial_cut(res,ges,r,r2,z,z2)
        print('Real after fiducial cut',len(res))
        print('GAN  after fiducial cut',len(ges))
        #resfit_r,resfit_g,er,eg = fit_res(real_amp,gan_amp,[4000,6500,500,3500], [5800,7400],[4000,6500,500,3500],[5800,7400],True,FileName,unit='A.U.')
        #real_res,gan_res,er,eg = cali_res(res,ges,[4000,6500,500,3500],[5000,7400],[4000,7800,500,3500],[6500,9000],FileName)
        ###real_res,gan_res,er,eg = cali_res(res,ges,[4000,6500,500,3500], [5000,7400],[4000,6500,500,3500],[5800,7400],FileName)
        '''
        plt.hist(r[:,0],bins=np.linspace(-200,200,40),label='X',alpha=0.5,histtype='step')
        plt.hist(r[:,1],bins=np.linspace(-200,200,40),label='Y',alpha=0.5,histtype='step')
        plt.hist(z,bins=np.linspace(-200,200,40),label='Z',alpha=0.5,histtype='step')
        plt.legend()
        plt.show()
        '''
        # quatratic correction
        #calibration_with_quadratic_fit(res,ges,FileName)
        with open("./pkls/gan_real_cali_E.pkl","rb") as file:
                cges = pickle.load(file)
                cres = pickle.load(file)
        #plot(gec,cges,rec,cres,2000,3500)
        #fit_res(cres,cges,[4000,2615,50,3500],[2000,2800],[4000,2615,50,3500],[2000,2800],True,FileName)
        #plot_corr(gec,cges,rec,cres)
        #rot_res(rec,cres,gec,cges)
        ### Here plot the rotated energy 
        #theta_min = 1.36 # without position correction
        theta_min = 1.193# GAN
        theta_min_r = 1.060 # Real
        fit_res(rec*np.sin(theta_min)+cres*np.cos(theta_min),gec*np.sin(theta_min)+cges*np.cos(theta_min),
            [4000,2615*np.sin(theta_min_r)+2615*np.cos(theta_min_r),50,3500],[2200*np.sin(theta_min_r)+2200*np.cos(theta_min_r),2800*np.sin(theta_min_r)+2800*np.cos(theta_min_r)],
            [4000,2615*np.sin(theta_min)+2615*np.cos(theta_min),50,3500],[2200*np.sin(theta_min)+2200*np.cos(theta_min),2800*np.sin(theta_min)+2800*np.cos(theta_min)],
            True,FileName+'-rotated','A.U.')
        
        
        '''
        print abs(eg - 2615)
        while abs(eg - 2615) > 0.001 * 2615:
            with open("./pkls/gan_real_cali_E.pkl","rb") as file:
                cges = pickle.load(file)
                cres = pickle.load(file)
            ## For phase 1
            #fit_res(cres,cges,[4000,2615,200,100],[1900,3000],[4000,2615,50,100],[1900,3500],True)
            ## For phase 2
            real_res,gan_res,er,eg = cali_res(cres,cges,[4000,2615,200,100],[1900,3000],[4000,2615,50,100],[1900,3000],FileName)
        '''
    else:
        import ROOT
        #RunNumber = glob.glob('/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase1/recon_*_gan_%s.root'%FileName)
        ROOT.gSystem.Load("libEXOUtilities")
        ROOT.gSystem.Load("libEXOFitting")
        gec,ges,rec,res = read(RunNumber,FileName,Phase)