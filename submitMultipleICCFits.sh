# MANTIDPYPATH="/home/ntv/mantid/mantid/bin"
MANTIDPYPATH="/SNS/users/ntv/mantid/mantid/release/bin"
NUMRUNS=5 #total number of runs



NUMOFFSET=1
NUMRUNSOFFSET=$(expr $NUMRUNS - $NUMOFFSET)

for i in $(eval echo "{0..$NUMRUNSOFFSET}")
do
    echo "running job with python doICCFit.py -r $i -p $MANTIDPYPATH"
    nohup python doICCFit.py -r $i -p $MANTIDPYPATH >/dev/null 2>&1  & 
done

