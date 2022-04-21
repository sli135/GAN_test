import os,glob
###################################################################
#
# phase1 lightmap: vanilla-base64
# phase2 lightmap: 2017_0nu-base64_v2
#
###################################################################
phase = 2
wf_filename = 'wgan_3d'#'critic-test'#'MPD-NoPooling-phase2-flatten'#'MPD-NoPooling-phase2-flatten_copy'#MPD-NoPooling-phase2-flatten-importAmpPenalty'#'MPD-NoPooling'#
wf_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_root/phase%i/%s/'%(phase,wf_filename)
wf = glob.glob(wf_path+'temp_full_*_gan_%s.root'%wf_filename)
#real = glob.glob(wf_path+'GAN_generated_waveform_*_real.csv'%(wf_filename))
print wf_path
#          use input rec wiregain v-wiregain cluster apdsignalfinder gridcorr lifecalib apdposcorr weighted-apd toutput
body ="""  use input rec wiregain v-wiregain cluster gridcorr lifecalib apdposcorr toutput
  /input/file %s
  /rec/LowerFitBoundWire 40
  /rec/UpperFitBoundWire 140
  /rec/ElectronicsDBFlavor timevartau
  /rec/enable_stage multiple_sig_finder true
  /rec/matched_filter_finder/APDSumThresholdFactor 5
  /rec/matched_filter_finder/APDSmoothWindow 4
  /rec/LowerFitBoundWireRestr 20
  /rec/UpperFitBoundWireRestr 30
  /rec/matched_filter_finder/DivideNoise false
  /rec/matched_filter_finder/UserVWireThreshold -1.0
  /rec/matched_filter_finder/VWireThresholdFactor 5.0
  /rec/collection_drift_velocity_mm_per_ns 0.00174
  /rec/drift_velocity_mm_per_ns 0.00171
  /rec/show_stage_status
  
  /rec/u_and_apd_fitter/verbose/outputtext apd false
  /rec/u_and_apd_fitter/verbose/outputtext u-wire false
  /rec/u_and_apd_fitter/verbose/outputtext v-wire false
  /rec/u_and_apd_fitter/verbose/outputtext induction false
  
  /rec/u_and_apd_fitter/verbose/plot apd false
  /rec/u_and_apd_fitter/verbose/plot u-wire false
  /rec/u_and_apd_fitter/verbose/plot v-wire false
  /rec/u_and_apd_fitter/verbose/plot induction false

  /rec/u_and_apd_fitter/verbose/plottofile apd false
  /rec/u_and_apd_fitter/verbose/plottofile u-wire false
  /rec/u_and_apd_fitter/verbose/plottofile v-wire false
  /rec/u_and_apd_fitter/verbose/plottofile induction false

  /wiregain/UWireDBFlavor pp_sourcecal_ci_timevar_2a
  /v-wiregain/VWireDBFlavor vanilla
  /cluster/drift_velocity_mm_per_ns 0.00171
  /cluster/ignore_induction True
  /cluster/collection_time 2940
  /gridcorr/gridcorrDBFlavor linear_expcorrections
  /lifecalib/manual 4500
  /apdposcorr/file /nfs/slac/g/exo/shaolei/GAN-test/GAN-test-lightmap/root_files/gan_lightmap_MPD-NoPooling-phase2-flatten.root
  /apdposcorr/skipCorr true
  /toutput/writeSignals false
  /toutput/onlyUWires false
  /toutput/trimWFs false
  /toutput/onlyDNNEvents false
  
  maxevents 100000
  printmodulo 1
  /toutput/file %s
  begin
"""
'''
/apdposcorr/file /nfs/slac/g/exo/shaolei/GAN-test/GAN-test-lightmap/root_files/gan_lightmap_MPD-phase2-flatten-AmpPenalty3.root
  /apdposcorr/skipCorr true
'''
recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/%s/'%(phase,wf_filename)
#recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/%s-lightmap-MPD-NoPooling-phase2-flatten-random/'%(phase,wf_filename)
#recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/%s/over2000/'%(phase,wf_filename)
#recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/%s/before-moving-baseline/'%(phase,wf_filename)
#recon_path = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/gan_lightmap_recon/MPD-phase2-flatten-AmpPenalty3/'%(phase)
try:
    os.mkdir(recon_path)
except:
    os.system('rm '+recon_path+'*')
for row in wf:
    realfile = row.replace('gan_%s'%wf_filename,'real')
    run = realfile[-18:-10]
    #outfile = '/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_recon/phase%i/gan_lightmap_recon/%s/recon_%s_gan_%s.root'%(phase,wf_filename,run,wf_filename) # gan_lightmap_recon/ or recon_with_lightmap/lightmap_2017_0nu-base64_v2
    outfile = recon_path+'recon_%s_gan_%s.root'%(run,wf_filename) 
    outfile_real = outfile.replace('gan_%s'%wf_filename,'real')

    filename = './output/gan_'+run+'.exo'
    f = open(filename,'w')
    f.write(body%(row,outfile))
    f.close()
    cmd = 'chmod 755 %s'%filename
    os.system(cmd)
    cmd = 'bsub -R rhel60 -W 72:00 -o output/gan_%s.out EXOAnalysis %s'%(run,filename)
    print cmd
    os.system(cmd)

    filename = './output/real_'+run+'.exo'
    f = open(filename,'w')
    f.write(body%(realfile,outfile_real))
    f.close()
    cmd = 'chmod 755 %s'%filename
    os.system(cmd)
    cmd = 'bsub -R rhel60 -W 72:00 -o output/real_%s.out EXOAnalysis %s'%(run,filename)
    print cmd
    os.system(cmd)


