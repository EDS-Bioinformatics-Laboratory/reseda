ping www.xs4all.nl

./mount-beehub.sh

cd git/progress/
git pull origin master
screen ./run.sh 

cd ../tbcell-miseq-pipeline/
git branch devel
git checkout devel
git config --global user.email "git@van-schaik.org"
git config --global user.name barbera
git pull origin devel

cp ~/SAMPLES-run07-* SAMPLES 
mv reference/* .
nano execute-all.sh 

nohup ./execute-all.sh > nohup.out 2> nohup.err < /dev/null &
