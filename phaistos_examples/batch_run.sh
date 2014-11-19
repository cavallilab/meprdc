#
# example shell script to perform multiple estimation sweeps using Phaistos and the inference module
# best run this from the build directory of phaistos, or make sure your paths are properly set. 
# 
#
BASE_PATH="YOUR PATH HERE"
SCRIPTS_PATH="PATH TO SCRIPTS"
NUM_THREADS=8
n=1
for j in {0..31};
do
  mkdir $BASE_PATH"/"$j"/"
  python $SCRIPTS_PATH"/makeconfig.py" $j $n $BASE_PATH"/run"$((j-1)).config
  ./phaistos --config-file $BASE_PATH"/run"$((j-1)).config  --threads $NUM_THREADS --energy-profasi-cached 1
  python  $SCRIPTS_PATH"/rdc_lambdas.py" $j `ls $BASE_PATH/$j/*_$n.dat`
done

