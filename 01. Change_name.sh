#!/bin/bash

## Author: Alicia Garcia-Roldan
## Date: October 2023

# This script is designed to rename the bins in a CheckM file. It is perfect for SqueezeMeta (v1.6.2post3) table 17.*.checkM located in intermediate folder.
# Now you can concatenate all the files you want without having repeated names for different bins. 
# You will need a files.txt (with the files you want to change) and samples.txt (with the names of the samples or the names you want to put in your new bins)

FILE=$1
LINENUMBER=1

#echo $SAMPLE

for FILE in $(cat files.txt)
do
	echo -e "Come on, let's live like cavemen [...] José Miguel, paint a buffalo on the wall..."
	
	SAMPLE=$(sed -n "${LINENUMBER}p" samples.txt)	
    	
	# The command works like: sed -n 1p samples.txt --> it is reading the first line of samples.txt

	echo $SAMPLE

	sed "s/^  concoct/${SAMPLE}_concoct/" ${FILE}  > ${SAMPLE}_var.checkM
	sed "s/^  maxbin/${SAMPLE}_maxbin/" ${SAMPLE}_var.checkM > ${SAMPLE}_var2.checkM
	sed "s/^  metabat2/${SAMPLE}_metabat2/" ${SAMPLE}_var2.checkM > ${SAMPLE}_var3.checkM

	rm ${SAMPLE}_var.checkM
	rm ${SAMPLE}_var2.checkM	
	mv ${SAMPLE}_var3.checkM ${SAMPLE}_var.checkM

	LINENUMBER=$((LINENUMBER + 1))	
done

# If you want to concatenate the files obtained

cat *_var.checkM > allbins_var.checkM

