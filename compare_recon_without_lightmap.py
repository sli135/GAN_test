import numpy as np 
import pickle
import csv

def read(phase):
    import ROOT
    ROOT.gSystem.Load("libEXOUtilities")
    ROOT.gSystem.Load("libEXOFitting")

    runNumbers = []
    with open('phase%i_runset.csv'%phase,'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if (row[1] == 'wTh') and row[3] == 'test' and row[0] > 3717:
                #print row[2]
                runNumbers.append(row[0])

    t = ROOT.TChain('tree','tree')
    for runNumber in runNumbers:
        print('Add file: %s'%runNumber)
        t.Add("/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/recon_with_lightmap/recon_%s*_gan_MPD-NoPooling.root" % (phase,runNumber)) # With lightmap
    ed = ROOT.EXOEventData()
    t.SetBranchAddress('EventBranch',ed)

    t2 = ROOT.TChain('tree','tree')
    for runNumber in runNumbers:
        print('Add file: %s'%runNumber)
        t2.Add("/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/recon_%s*_gan_MPD-NoPooling.root" % (phase,runNumber))
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

    for i in range(t.GetEntries()):
        t.GetEntry(i)
        t2.GetEntry(i)
        nev = ed.fEventNumber
        nev2 = ed2.fEventNumber

        ncc = ed.GetNumChargeClusters()
        ncc2 = ed2.GetNumChargeClusters()
        if ncc != ncc2:
            wrong_event_2.append([i,nev,nev2,ncc,ncc2])
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

        nsc = ed.GetNumScintillationClusters()
        nsc2 = ed2.GetNumScintillationClusters()
        if nsc != nsc2:
            wrong_event.append([i,nev,nev2,nsc,nsc2])
        esc = 0
        for j in range(nsc):
            sc = ed.GetScintillationCluster(j)
            esc += sc.fRawEnergy
        gan_e_scint.append(esc)
        esc = 0
        for j in range(nsc2):
            sc = ed2.GetScintillationCluster(j)
            esc += sc.fRawEnergy
        if esc == 6756.9145724690461 and filename == 'smallPhase1se':
            continue
        real_e_scint.append(esc)

        Id.append(i)
        Id2.append(i)

    with open("./pkls/recon_compare_lightmap_E.pkl","wb") as file:
        pickle.dump(np.array(gan_e_scint),file)
        pickle.dump(np.array(real_e_scint),file)

def main():
    import matplotlib.pyplot as plt
    with open("./pkls/recon_compare_lightmap_E.pkl","rb") as file:
        lightE_with_lightmap = pickle.load(file)
        lightR_without_lightmap = pickle.load(file)
            
    left = 500
    right = 12000
    nbins = (right - left) / 40
    mybins = np.linspace(left,right,nbins+1)
    e_step = (right-left)/nbins
    nGANeScint, binsGANeScint = np.histogram(lightE_with_lightmap,bins=mybins)
    nRealeScint, binsRealeScint = np.histogram (lightR_without_lightmap,bins=mybins)

    err_ges = get_errs(nGANeScint)
    err_res = get_errs(nRealeScint)
    cent  = get_centers(mybins,e_step)
    plt.rcParams['figure.figsize'] = [12, 10]
    plt.errorbar(cent,nGANeScint,yerr=err_ges,ecolor='lightgrey',fmt='o',color='r',markersize=6.,label='light E with lightmap')
    plt.errorbar(cent,nRealeScint,yerr=err_res,ecolor='lightgrey',fmt='o',color='b',markersize=6.,label='light E without lightmap')
    plt.legend()
    plt.title(r'light E',fontsize=20)
    plt.xlabel(r'fRawEnergy [A.U.]',fontsize=20)
    plt.ylabel(r'Counts/(%.0f [A.U.])'%(e_step),fontsize=20)
    plt.savefig('pics/light_E_with_light_gan.png',bbox_inches = 'tight')
    plt.show()
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
if __name__ == '__main__':
    #read(1)
    main()