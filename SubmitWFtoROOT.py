import os,glob

phase = 2
wf_filename = 'wgan_3d'#'critic-test'#'MPD-NoPooling-phase2-flatten_copy'#'MPD-NoPooling-phase2-flatten'#'phase1_se_reweight'
wf_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_wf/gan/phase%i/'%phase
root_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_root/phase%i/%s/'%(phase,wf_filename)
try:
    os.mkdir(root_path)
except:
    os.system('rm '+root_path+'*')

wf = glob.glob(wf_path+'GAN_generated_waveform_*_%s.csv'%wf_filename)
body = 'bsub -R rhel60 -W 72:00 -o output/make_wf_%s.out python myGAN_wf_root.py --Phase %i --RunNumber %s --FileNumber %s --WFfile %s'
i = -len(wf_filename) - 6
for row in wf:
    runNumber,fileNumber = row[i - 7:i - 3],row[i]
    cmd = body % (runNumber,phase,runNumber,fileNumber,wf_filename)
    print cmd
    os.system(cmd)