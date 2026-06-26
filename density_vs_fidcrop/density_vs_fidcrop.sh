echo "Run ${0} for ${1} and ncrops=${2}\n"

LOCAL=.
#def analysis prgm
PRGM=$LOCAL/density_vs_fidcrop.py

#def target dir
SIMUS=$LOCAL/simus_sample

#def save directory
SAVE=$LOCAL/pixs
# SAVE=$LOCAL/pixs/fid_energy_script_zoom_20_keV

#check analysis prgm
if [ ! -f $PRGM ]
then
  echo "no analysis prgm. Script abort"
  exit 1
fi

#check target dir
if [ ! -d $SIMUS ]
then
  echo "no target directory. Script abort"
  exit 1
fi

#check save dir
if [ ! -d $SAVE ]
then
  echo "no save directory. Script abort"
  exit 1
fi

################################################################################
echo "Matching files found : "
for file in $SIMUS/*${1}.root
do
  if [ -f $file ]
  then
    echo "|\t $file"
  else
    echo "found no file such as $SIMUS/*${1}.root.\n\nEnd of prgm"
    exit 1
  fi
done

#loop over the file in $SIMUS
for file in $SIMUS/*${1}.root
do
  if [ -f $file ]
  then
    echo '--------------------------------------------------------------------------------'
    python3 $PRGM $file $SAVE ${2}
    echo '--------------------------------------------------------------------------------'
  fi
done
