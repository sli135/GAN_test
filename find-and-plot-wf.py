import numpy as np 
import os 
import argparse
import glob


def show_spectrum(runNumber,filename,phase):
    
    recon_gan_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/%s/recon_%s_gan_%s.root'%(phase,filename,runNumber,filename)
    recon_real_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/%s/recon_%s_real.root'%(phase,filename,runNumber)
    t_gan,t_real = ROOT.TChain('tree','tree'),ROOT.TChain('tree','tree')
    t_gan.Add(recon_gan_path)
    t_real.Add(recon_real_path)
    ed_gan,ed_real = ROOT.EXOEventData(),ROOT.EXOEventData()
    t_gan.SetBranchAddress('EventBranch',ed_gan)
    t_real.SetBranchAddress('EventBranch',ed_real)
    print('GAN tree',t_gan.GetEntries())
    print('Real tree',t_real.GetEntries())
    sample_targets = np.array(range(t_gan.GetEntries()))
    np.random.shuffle(sample_targets)
    save_path = './pkls/wf_samples/low-test-sample-%s-*.pkl'%(filename)
    print save_path
    os.system('rm '+save_path)
    
    for i in sample_targets[:100]:
        t_gan.GetEntry(i)
        t_real.GetEntry(i)
        nsc_gan,nsc_real = ed_gan.GetNumScintillationClusters(),ed_real.GetNumScintillationClusters()
        esc_gan,esc_real = 0,0 
        for j in range(nsc_gan):
            sc = ed_gan.GetScintillationCluster(j)
            esc_gan += sc.fRawEnergy
        for k in range(nsc_real):
            sc = ed_real.GetScintillationCluster(k)
            esc_real += sc.fRawEnergy
        body = 'bsub -R rhel60 -W 72:00 -o output/wf_samples/%s-%i.out python low-test-plot.py --RunNumber %s --nev %i  --FileName %s --Phase %i'
        cmd = body%(filename,i,runNumber,i,filename,phase)
        print(cmd)
        os.system(cmd)
def plot_all_wf(runNumber,filename,phase):
    pickle_path = '/nfs/slac/g/exo/shaolei/GAN-test/pkls/wf_samples/'
    files = glob.glob(pickle_path+'low-test-sample-'+filename+'-*.pkl')
    start = len(pickle_path+'low-test-sample-'+filename) + 1
    samples = [f[start:-4]for f in files]
    for i in samples:
        body = 'bsub -R rhel60 -W 72:00 -o output/wf_samples/%s-%s.out python low-test-plot.py --RunNumber %s --nev %s  --FileName %s --Phase %i --ReadFile t'
        cmd = body%(filename,i,runNumber,i,filename,phase)
        print(cmd)
        os.system(cmd)
if __name__ == '__main__':
    import ROOT
    ROOT.gSystem.Load("libEXOUtilities")
    ROOT.gSystem.Load("libEXOFitting")
    # Phase 2
    name = 'critic-test'#'MPD-NoPooling-phase2-flatten_copy'
    show_spectrum('8673-000',name,2)
    #plot_all_wf('8673-000',name,2)
    # Phase 1
    #show_spectrum('8673-000','MPD-phase2-flatten-AmpPenalty5',1)
    #plot_all_wf('8673-000','MPD-phase2-flatten-AmpPenalty5',1)
