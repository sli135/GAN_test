##############################################################
#
# python PlotDeadAPD.py --RunNumber 5872 --ModelName scintE
#
##############################################################
import numpy as np
import argparse
import matplotlib as mpl
mpl.use('Agg')
#mpl.rcParams['mathtext.fontset'] = 'cm'
#mpl.rc('font', family='serif')

import matplotlib.pyplot as plt 
import matplotlib.image as img 
from matplotlib import ticker, cm

import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator

plt.rcParams['figure.figsize'] = [10, 5]

import sys, code, os, glob, csv
import math
import random

from scipy import stats
from scipy import interpolate 
from scipy.optimize import curve_fit

from array import array

import pickle

import time
def getEs(apd_waves):
    totalE = 0
    EperCh1 = []
    EperCh2 = []
    bmeans = []
    brmss = []
    if len(apd_waves) != 74:
        print('too bad!')
        return -1, EperCh
    for i in range(0,74):
        bmean = np.mean(apd_waves[i][0:50])
        brms = np.std(apd_waves[i][0:50])
        bmeans.append(bmean)
        brmss.append(brms)
        e_i = np.max(apd_waves[i])
        if e_i < 3.*brms:
            e_i = 0 
        if i<37:
            EperCh1.append(e_i)
        else:
            EperCh2.append(e_i)
        totalE = totalE+e_i
    return totalE, np.asarray(EperCh1), np.asarray(EperCh2), bmeans, brmss  

xAPDcenter = 0 
yAPDcenter = 0
nAPDrows = 21 #number of rows of APDs
APDrowsLeft  = np.array([2,3,5,5,7,7,9,9,9,8,9,8,8,8,8,7,8,7,5,2,0]) # APDs in each row
APDrowsRight = np.array([0,2,4,7,7,7,7,8,7,8,8,8,8,9,8,7,6,5,4,3,2]) # APDs in each row    
APDnumber = 0

APD_SPACING = 2.218182 #cm
APDRADIUS = 0.8 #cm
NUMBER_APD_PER_PLANE=259
NUMBER_APD_CHANNELS_PER_PLANE = 37

fAPDx = np.zeros(NUMBER_APD_PER_PLANE)
fAPDy = np.zeros(NUMBER_APD_PER_PLANE)


for i in range(0,nAPDrows):
    for j in range(0, (APDrowsLeft[i]+1)):
        if APDrowsRight[i] != 0:
            for k in range(0, APDrowsRight[i]):
                if (j != APDrowsLeft[i]) and (k == 0):
                    fAPDy[APDnumber] = ((nAPDrows + 1.)/2. - (i + 1.)) * APD_SPACING*np.sqrt(3.)/2.
                    if APDrowsLeft[i] == 2 and APDrowsRight[i] == 0:
                        fAPDx[APDnumber] = (-APDrowsLeft[i] + j) * APD_SPACING
                    else:
                        if np.fmod(i,2.) > 0.5:
                            fAPDx[APDnumber] = (-(APDrowsLeft[i])+ 0.5 + j) * APD_SPACING
                        else:
                            fAPDx[APDnumber] = (-(APDrowsLeft[i]) + j + 1.) * APD_SPACING
                    APDnumber = APDnumber+1
                elif j == APDrowsLeft[i]:
                    fAPDy[APDnumber] = ((nAPDrows + 1.)/2. - (i + 1.)) * APD_SPACING*np.sqrt(3.)/2.
                    if np.fmod(i,2.) > 0.5:
                        fAPDx[APDnumber] = (0.5 + k) * APD_SPACING
                    else:
                        fAPDx[APDnumber] = (1. + k) * APD_SPACING
                    APDnumber = APDnumber + 1
        else:
            if j != APDrowsLeft[i]:
                fAPDy[APDnumber] = ((nAPDrows + 1.)/2. - (i + 1.)) * APD_SPACING*np.sqrt(3.)/2.
                if APDrowsLeft[i] == 2 and APDrowsRight[i] == 0:
                    fAPDx[APDnumber] = (-APDrowsLeft[i] + j) * APD_SPACING
                else:
                    if np.fmod(i,2.) > 0.5:
                        fAPDx[APDnumber] = (-(APDrowsLeft[i])+ 0.5 + j) * APD_SPACING
                    else:
                        fAPDx[APDnumber] = (-(APDrowsLeft[i]) + j + 1.) * APD_SPACING
                APDnumber = APDnumber + 1
      
    
print(APDnumber)

theta = (12-30.)*np.pi/180.
cost = np.cos(theta)
sint = np.sin(theta)
fAPDxR = np.zeros(NUMBER_APD_PER_PLANE)
fAPDyR = np.zeros(NUMBER_APD_PER_PLANE)

for i in range(0,APDnumber):
    x0 = fAPDx[i]
    y0 = fAPDy[i]
    x  = cost*x0 + sint*y0
    y  = -1.*sint*x0 + cost*y0
    fAPDxR[i]=x
    fAPDyR[i]=y
APDPlane1_gang_table = np.array([1,1,
                    1,1,1,2,2,
                    34,34,1,1,2,2,2,3,3,
                    34,34,34,4,4,2,2,3,3,3,7,7,
                    33,33,34,34,4,4,4,5,5,3,3,7,7,7,
                    33,33,33,36,36,4,4,5,5,5,10,10,7,7,
                    32,32,33,33,36,36,36,6,6,5,5,10,10,10,8,8,
                    32,32,32,35,35,36,36,6,6,6,12,12,10,10,8,8,8,
                    32,32,35,35,35,37,37,6,6,12,12,12,11,11,8,8,
                    27,27,35,35,37,37,37,31,31,12,12,11,11,11,9,9,
                    27,27,27,29,29,37,37,31,31,31,18,18,11,11,9,9,9,
                    27,27,29,29,29,30,30,31,31,18,18,18,16,16,9,9,
                    26,26,29,29,30,30,30,24,24,18,18,16,16,16,13,13,
                    26,26,26,28,28,30,30,24,24,24,17,17,16,16,13,13,13,
                    26,26,28,28,28,23,23,24,24,17,17,17,14,14,13,13,
                    25,25,28,28,23,23,23,22,22,17,17,14,14,14,
                    25,25,25,21,21,23,23,22,22,22,15,15,14,14,
                    25,25,21,21,21,20,20,22,22,15,15,15,
                    21,21,20,20,20,19,19,15,15,
                    20,20,19,19,19,
                    19,19])
  

APDPlane2_gang_table = np.array([20,20,
                    20,20,20,21,21,
                    16,16,20,20,21,21,21,22,22,
                    16,16,16,23,23,21,21,22,22,22,26,26,
                    15,15,16,16,23,23,23,24,24,22,22,26,26,26,
                    15,15,15,18,18,23,23,24,24,24,29,29,26,26,
                    14,14,15,15,18,18,18,25,25,24,24,29,29,29,27,27,
                    14,14,14,17,17,18,18,25,25,25,31,31,29,29,27,27,27,
                    14,14,17,17,17,19,19,25,25,31,31,31,30,30,27,27,
                    9,9,17,17,19,19,19,13,13,31,31,30,30,30,28,28,
                    9,9,9,11,11,19,19,13,13,13,37,37,30,30,28,28,28,
                    9,9,11,11,11,12,12,13,13,37,37,37,35,35,28,28,
                    8,8,11,11,12,12,12,6,6,37,37,35,35,35,32,32,
                    8,8,8,10,10,12,12,6,6,6,36,36,35,35,32,32,32,
                    8,8,10,10,10,5,5,6,6,36,36,36,33,33,32,32,
                    7,7,10,10,5,5,5,4,4,36,36,33,33,33,
                    7,7,7,3,3,5,5,4,4,4,34,34,33,33,
                    7,7,3,3,3,2,2,4,4,34,34,34,
                    3,3,2,2,2,1,1,34,34,
                    2,2,1,1,1,
                    1,1])

def get_gang_number(apd,plane):
    if plane == 0:
        gang = APDPlane1_gang_table[apd] - 1
    elif plane == 1:
        gang = APDPlane2_gang_table[apd] - 1
    else:
        gang = -999
    return gang

def get_apds_from_gang(gang_no,plane):
    apds = []
    for i in range(0,APDnumber):
        g_i = get_gang_number(i,plane)
        if g_i == gang_no:
            apds.append(i)
    return apds
def get_color(EperCh,i,enorm=-1):
    if enorm==-1:
        enorm = np.max(EperCh)
    shade = (enorm-EperCh[i])/enorm
    return (shade,shade,shade)
def plot_one_event(event,savename=''):

    plt.rcParams['figure.figsize'] = [10,5]
 
    fig1 = plt.figure()
    frame1=fig1.add_axes([0.,0.,0.5,1])
    frame1.set_xlabel('X (cm)',fontweight='bold',fontsize=20)
    frame1.set_ylabel('Y (cm)',fontweight='bold',fontsize=20)
    frame1.axis((-25,25,-25,25))
    frame1.tick_params(labelsize=14)
    
    frame2=fig1.add_axes([0.5,0,0.5,1])
    frame2.set_xlabel('X (cm)',fontweight='bold',fontsize=20)
    #frame2.set_ylabel('Y (cm)',fontweight='bold',fontsize=20)
    frame2.axis((-25,25,-25,25))
    frame2.tick_params(labelsize=14)

    totalE, EperCh1, EperCh2,bmeans, brmss = getEs(event)
    
    for i in range(0,len(EperCh1)): 
        apds1 = get_apds_from_gang(i,0)
        apds2 = get_apds_from_gang(i,1)
        for j in range(0,len(apds1)):
            circ = plt.Circle((fAPDxR[apds1[j]],fAPDyR[apds1[j]]),APDRADIUS,
                              color=get_color(EperCh1,i),fill=True,lw=0.5,ls='--')
            frame1.add_artist(circ)     
        for j in range(0,len(apds2)):
            circ = plt.Circle((fAPDxR[apds2[j]],fAPDyR[apds2[j]]),APDRADIUS,
                              color=get_color(EperCh2,i),fill=True,lw=0.5,ls='--')
            frame2.add_artist(circ)
    
    frame1.text(.05, .95, "Plate 1 (S8)", transform=frame1.transAxes, ha="left", va="top")       
    frame2.text(.05, .95, "Plate 2 (S2)", transform=frame2.transAxes, ha="left", va="top")   
    
    if savename!='':
        plt.savefig(mypath+'./'+savename+'.png',bbox_inches = 'tight')
def plot_all_events(waves,cutoff=10,savename=''):
    plt.rcParams['figure.figsize'] = [10,5]
 
    fig1 = plt.figure()
    frame1=fig1.add_axes([0.,0.,0.5,1])
    frame1.set_xlabel('X (cm)',fontweight='bold',fontsize=20)
    frame1.set_ylabel('Y (cm)',fontweight='bold',fontsize=20)
    frame1.axis((-25,25,-25,25))
    frame1.tick_params(labelsize=14)
    
    frame2=fig1.add_axes([0.5,0,0.5,1])
    frame2.set_xlabel('X (cm)',fontweight='bold',fontsize=20)
    #frame2.set_ylabel('Y (cm)',fontweight='bold',fontsize=20)
    frame2.axis((-25,25,-25,25))
    frame2.tick_params(labelsize=14)
    
    EperCh1_tot = np.zeros(NUMBER_APD_CHANNELS_PER_PLANE)
    EperCh2_tot = np.zeros(NUMBER_APD_CHANNELS_PER_PLANE)
    for k in range(0,cutoff):
        totalE, EperCh1, EperCh2,bmeans, brmss = getEs(waves[k])
        for w in range(0,len(EperCh1)):
            EperCh1_tot[w] = EperCh1_tot[w]+EperCh1[w]
            EperCh2_tot[w] = EperCh2_tot[w]+EperCh2[w]
    enorm = np.max(EperCh1_tot)
    if np.max(EperCh2_tot) > enorm:
        enorm = np.max(EperCh2_tot)

    shadow1 = [(get_color(EperCh1_tot,i,enorm)[0],EperCh1_tot[i],i) for i in range(36)]
    shadow1 = sorted(shadow1, key = lambda x: x[1])
    shadow2 = [(get_color(EperCh2_tot,i,enorm)[0],EperCh2_tot[i],i) for i in range(36)]
    shadow2 = sorted(shadow2, key = lambda x: x[1])
    print(shadow1)
    print(shadow2)
    
    for i in range(0,len(EperCh1)): 
        apds1 = get_apds_from_gang(i,0)
        apds2 = get_apds_from_gang(i,1)
        for j in range(0,len(apds1)):
            circ = plt.Circle((fAPDxR[apds1[j]],fAPDyR[apds1[j]]),APDRADIUS,
                              color=get_color(EperCh1_tot,i,enorm),fill=True,lw=0.5)
            frame1.add_artist(circ)
        for j in range(0,len(apds2)):
            circ = plt.Circle((fAPDxR[apds2[j]],fAPDyR[apds2[j]]),APDRADIUS,
                              color=get_color(EperCh2_tot,i,enorm),fill=True,lw=0.5)
            frame2.add_artist(circ)
            
    frame1.text(.05, .95, "Plate 1 (S8)", transform=frame1.transAxes, ha="left", va="top")       
    frame2.text(.05, .95, "Plate 2 (S2)", transform=frame2.transAxes, ha="left", va="top")  

    if savename!='':
        frame1.text(.45, .95, savename, transform=frame1.transAxes, ha="left", va="top")
        frame2.text(.45, .95, savename, transform=frame2.transAxes, ha="left", va="top")
        plt.savefig('./pics/dead-apds/'+savename+'.png',bbox_inches = 'tight')
        plt.clf()
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--RunNumber', type=int, nargs='+', default=6230)
    parser.add_argument('--ModelName', type=str, nargs=1, default='scintE')
    parser.add_argument('--Phase', type=int, default=1)
    args = parser.parse_args()
    runNumber = args.RunNumber[0]
    modelname = args.ModelName[0]
    phase = args.Phase
    mypath = "/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_wf"

    wf_gan = glob.glob(mypath+'/gan/phase%i/GAN_generated_waveform_%i-*_%s.csv'%(phase,runNumber,modelname))
    wf_real = glob.glob(mypath+'/real/phase%i/real_waveform_%i-*.csv'%(phase,runNumber))
    wfs = []
    for file in wf_gan:
        with open(file,'r') as f:
            for row in f:
                x = np.fromstring(row, dtype=float, sep=',')
                wfs.append(x)
    wfs = np.array(wfs)
    print(wfs.shape)
    wfs = wfs.reshape(len(wfs),74,350)
    print 'GAN events',wfs.shape
    plot_all_events(wfs,len(wfs),savename='gan_%i'%runNumber)
    wfs = []
    for file in wf_real:
        with open(file,'r') as f:
            for row in f:
                x = np.fromstring(row, dtype=float, sep=',')
                wfs.append(x)
    wfs = np.array(wfs)
    wfs = wfs.reshape(len(wfs),74,350)
    print 'Real events',wfs.shape
    plot_all_events(wfs,len(wfs),savename='real_%i'%runNumber)
    '''
    for chn in range(74):
        for i in wfs[:100]:
            plt.plot(i[chn])
        plt.title('chn %i' % chn)
        plt.show() # dead chn 11 26 39 53
    '''



