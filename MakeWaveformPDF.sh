pdfjam --a4paper --fitpaper true pics/wf_samples/waveforms_test_*.png
rm -rf pics/wf_samples/waveforms_test_*.png
mv waveforms_test_*pdf pics/wf_samples/ 

pdfjam --a4paper --fitpaper true pics/wf_samples/waveforms_platter_*.png
rm -rf pics/wf_samples/waveforms_platter_*.png
mv waveforms_platter_*pdf pics/wf_samples/ 