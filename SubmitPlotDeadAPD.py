import os,glob

phase = 2
wf_filename = 'wgan_3d'#'MPD-phase2-flatten-AmpPenalty3'#'phase1_se_reweight'
wf_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_wf/gan/phase%i/'%phase
wf = glob.glob(wf_path+'GAN_generated_waveform_*_%s.csv'%wf_filename)
body = 'bsub -R rhel60 -W 72:00 -o output/apd_%s.out python PlotDeadAPD.py --RunNumber %s --ModelName %s --Phase %i'
i = -len(wf_filename) - 6
for row in wf:
    runNumber,fileNumber = row[i - 7:i - 3],row[i]
    cmd = body % (runNumber,runNumber,wf_filename,phase)
    print cmd
    os.system(cmd)