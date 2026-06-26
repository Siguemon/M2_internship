
LOCAL=.

PRGM=$LOCAL/main.py

SIMULATIONS=$LOCAL/data/simus #the path for the simulations
SIMU_NAMES=$LOCAL/data/simu_names.json #the path for the file that lists the simulations (the one we will use to loop on)
OUT_FILEDIR=$LOCAL/data/out #the path on which the output will go

python3 $PRGM $SIMULATIONS $SIMU_NAMES $OUT_FILEDIR/$1
