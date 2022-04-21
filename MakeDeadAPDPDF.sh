pdfjam --a4paper --fitpaper true pics/dead-apds/gan_*.png
rm -rf pics/dead-apds/gan_*.png
mv gan_*pdf pics/dead-apds/ 
pdfjam --a4paper --fitpaper true pics/dead-apds/real_*.png
rm -rf pics/dead-apds/real_*.png
mv real_*pdf pics/dead-apds/ 