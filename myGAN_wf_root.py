#################################################################################
#
# python myGAN_wf_root.py --RunNumber 5872 --FileNumber 0 --WFfile [gan_wf/*.csv] unflatten
#
#################################################################################

import os, re, sys, glob, math, random, pickle
import numpy as np
import ROOT
import argparse

def merge_data(runNumber,wv_file,fileNumber,phase):
	# Some constants
	rotationAngle = 0.5312
	time_of_events = 1357754448
	nchannels = 226
	subtractBaseline = True
	#Create ROOT file and make the Event Data Branch
	standard_wf_length = 2048
	standard_period = 1*ROOT.CLHEP.microsecond
	print("*************************************************************")
	print("RunNumber:       %i"%runNumber)
	print("*************************************************************")
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
	preprocessedFiles = '/nfs/slac/g/exo_data3/exo_data/data/WIPP/masked/5872/masked00005872-000.root' #'./run_5872_tree.root'

	print('Opening preprocessed file %s'%(preprocessedFiles))
	fProcessed = ROOT.TFile(preprocessedFiles,'READ')
	tProcessed = fProcessed.Get('tree')#dataTree')
	print("tProcessed events:",tProcessed.GetEntries())

	print('Opening waveform files %s'%(waveformFiles))
	tWF = ROOT.TChain('tree','tree')
	tWF.Add(waveformFiles)
	print("tWF events:",tWF.GetEntries())
	'''
	#### import GAN wv #############
	wv_path = "/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_wf/gan/phase%i/"%phase
	try:
		wv_name = "GAN_generated_waveform_%i-%03i_%s.csv" % (runNumber,fileNumber,wv_file)
	except:
		wv_name = "GAN_generated_waveform_%i_%s.csv" % (runNumber,wv_file)
	xdata = []
	with open(wv_path + wv_name,"rb") as f:
		for line in f:
			x = np.fromstring(line, dtype=float, sep=',')
			xdata.append(x[0:74*350])
	f.close()
	xdata = np.array(xdata)
	xdata = xdata.reshape(xdata.shape[0],74,350)
	print("waveform No:",xdata.shape[0])
	#### Save root file ############
	my_save_path = "/nfs/slac/g/exo_data8/exo_data/users/shaolei/gan_root/phase%i/%s/" % (phase,wv_file)
	open_file_real = ROOT.TFile(my_save_path+"temp_full_%i-%03i_real.root" % (runNumber,fileNumber), "recreate")
	save_tree_real = ROOT.TTree("tree", "tree")
	ed_real = ROOT.EXOEventData()
	save_tree_real.Branch("EventBranch", ed_real)
	open_file_gan = ROOT.TFile(my_save_path+"temp_full_%i-%03i_gan_%s.root" % (runNumber,fileNumber,wv_file), "recreate")
	save_tree_gan = ROOT.TTree("tree", "tree")
	ed_gan = ROOT.EXOEventData()
	save_tree_gan.Branch("EventBranch", ed_gan)

	n = xdata.shape[0]
	ni = 0

	for i in range(0,tWF.GetEntries()):
		#print('i=%i, ni=%i'%(i,ni))
		if ni == n:
			break

		tProcessed.GetEntry(i)
		#rec_ED = tProcessed.EventBranch#Summary
		tWF.GetEntry(i)
		#raw_ED = tWF.EventBranch
		wf_data = raw_ED.GetWaveformData()

		#if rec_ED.energy<700. or rec_ED.isNoise or rec_ED.isVetoed :
		ncl = rec_ED.GetNumChargeClusters()
		eq = rec_ED.GetTotalPurityCorrectedEnergy()
		nsc = rec_ED.GetNumScintillationClusters()

		if i % 1000 == 0:
			print('%i events processed, %i event accepted'%(i,ni))

		if ncl!=1 or nsc!=1 or wf_data.fNumSamples!=2048: # or eq<2200. 
			continue
		cc = rec_ED.GetChargeCluster(0)
		sc = rec_ED.GetScintillationCluster(0)
		x,y,z,e_charge,e_scint = cc.fX,cc.fY,cc.fZ,cc.fPurityCorrectedEnergy,sc.fRawEnergy
		if x ** 2 + y ** 2 >= 162 ** 2 or abs(z) >= 180 or e_charge <= 600 or e_scint <= 0:
			continue


		wf_data.Decompress()
		num_waves  = wf_data.GetNumWaveforms()
		ed_real.fRunNumber = raw_ED.fRunNumber
		ed_real.fEventNumber = raw_ED.fEventNumber
		ed_real.fEventHeader = raw_ED.fEventHeader
		ed_gan.fRunNumber = raw_ED.fRunNumber
		ed_gan.fEventNumber = raw_ED.fEventNumber
		ed_gan.fEventHeader = raw_ED.fEventHeader
		minSample = 850
		t_scint = sc.fTime / 1000 - rec_ED.fEventHeader.fTriggerOffset
		if t_scint > -180:
			minSample = minSample + int(t_scint)
		for k in range(num_waves):
			#Create EXO Waveform from numpy array
			waveform = wf_data.GetWaveform(k)
			chan = waveform.fChannel
			arr_real = np.array(waveform)
			arr_real = arr_real.astype(np.double)
			arr_gan = np.array(waveform)
			arr_gan = arr_gan.astype(np.double)
			# dead chn 11 26 39 53
			if ROOT.EXOMiscUtil.TypeOfChannel(chan) == ROOT.EXOMiscUtil.kAPDGang: #and chan - 152 not in [11,26,39,53]:
				mu, sigma = np.mean(arr_gan[:50]), np.std(arr_gan[:50])
				mu2 = np.mean(xdata[ni][chan - 152][:40])
				data = xdata[ni][chan - 152] - mu2 + mu
				noise_1 = arr_gan[:minSample]
				noise_2 = arr_gan[minSample+350:]#arr_gan[:2048-minSample-350]
				arr_gan = np.concatenate([noise_1,data])
				arr_gan = np.concatenate([arr_gan,noise_2])
			dw_wf_real = ROOT.EXODoubleWaveform(arr_real, len(arr_real))
			dw_wf_real.SetSamplingPeriod(standard_period)
			exo_wf_real = ed_real.GetWaveformData().GetNewWaveform()
			exo_wf_real.ConvertFrom(dw_wf_real)
			exo_wf_real.fChannel = chan
			dw_wf_real.IsA().Destructor( dw_wf_real )
			dw_wf_gan = ROOT.EXODoubleWaveform(arr_gan, len(arr_gan))
			dw_wf_gan.SetSamplingPeriod(standard_period)
			exo_wf_gan = ed_gan.GetWaveformData().GetNewWaveform()
			exo_wf_gan.ConvertFrom(dw_wf_gan)
			exo_wf_gan.fChannel = chan
			dw_wf_gan.IsA().Destructor( dw_wf_gan )

		ed_real.GetWaveformData().fNumSamples = standard_wf_length
		ed_real.GetWaveformData().Compress()
		save_tree_real.Fill() #Fill tree with ED
		ed_real.Clear() #Clear ED for next event
		ed_gan.GetWaveformData().fNumSamples = standard_wf_length
		ed_gan.GetWaveformData().Compress()
		save_tree_gan.Fill() #Fill tree with ED
		ed_gan.Clear() #Clear ED for next event

		ni = ni+1

	if not save_tree_gan:
		print('Did not create tree successfully.')
	else:
		print('tree Entries:',save_tree_gan.GetEntries())
	open_file_real.cd()
	save_tree_real.Write() #Write to the file
	open_file_real.Close()
	open_file_real.IsA().Destructor( open_file_real )
	open_file_gan.cd()
	save_tree_gan.Write() #Write to the file
	open_file_gan.Close()
	open_file_gan.IsA().Destructor( open_file_gan )


if __name__ == "__main__":
	ROOT.gSystem.Load("libEXOUtilities")
	ROOT.gSystem.Load("libEXOFitting")
	parser = argparse.ArgumentParser()
	parser.add_argument('--RunNumber', type=int, nargs=1, default=5872)
	parser.add_argument('--FileNumber', type=int, nargs=1,default=0)
	parser.add_argument('--WFfile', type=str)
	parser.add_argument('--Phase', type=int,default = 1)

	args = parser.parse_args()
	RunNumber = args.RunNumber[0]
	wv_file = args.WFfile
	FileNumber = args.FileNumber[0]
	phase = args.Phase
	merge_data(RunNumber,wv_file,FileNumber,phase)
