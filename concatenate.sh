#
cp SAMPLES-basta SAMPLES ; wait
./copy-from-beehub.sh ; wait
myfiles=`cat LOCAL_SAMPLES`
python ConcatenateCloneFiles.py -r 20210401-DataSheet-BASTA.json -n BASTA -c IGH_HUMAN -pre vjcdr3-clones-mut- $myfiles ; wait
rm $myfiles
#
echo 'FINISHED'
