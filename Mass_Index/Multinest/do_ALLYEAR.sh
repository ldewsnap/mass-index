rdir="FULL_MOVIE"
rdir_ORB="FULL_MOVIE_ORB"
data_dir="data_FULL"
mkdir $rdir_ORB

mkdir $rdir
mkdir $data_dir


touch ./$data_dir/FULL_ORB.res
rm  ./$data_dir/FULL_ORB.res

touch ./$data_dir/FULL_ORB.res
rm  ./$data_dir/FULL_ORB.res



#for fil in ./database/QUA

for long in $(seq 0 1 0)
do
echo "" > quad_temp
echo "" >> quad_temp
echo "" >> quad_temp
echo "2011 2012 2013 2014 2015" >> quad_temp
echo "0" >> quad_temp
echo "360" >> quad_temp

cat quad_temp | sh picker.sh 
flnm=$long
echo $flnm
sh mass.sh
convert -quality 100 -density 300 -rotate 90 -flatten plot.eps ./$rdir/$flnm""_MEV.png
convert -quality 100 -density 300 -rotate 90 -flatten pgplot.ps ./$rdir/$flnm""_MEV_cov.png
cp datafile ./$data_dir/$long
awk '{ print '$long' "," $0 "," '`wc -l datafile | awk '''{print $1}'''`' }'  < results.out >> ./$data_dir/FULL.res



#cat quad_temp | sh picker_ORB.sh 
#c=`wc datafile -l | awk '{print $1}'`;  if [ "$c" -gt "10" ]; then
#flnm=$long
#echo $flnm
#sh mass.sh
#convert -quality 100 -density 300 -rotate 90 -flatten plot.eps ./$rdir_ORB/$flnm""_MEV.png
#convert -quality 100 -density 300 -rotate 90 -flatten pgplot.ps ./$rdir_ORB/$flnm""_MEV_cov.png
#cp datafile ./$data_dir/$long""_ORB
#awk '{ print '$long' "," $0 "," '`wc -l datafile | awk '''{print $1}'''`' }'  < results.out >> ./$data_dir/FULL_ORB.res
#fi

done
