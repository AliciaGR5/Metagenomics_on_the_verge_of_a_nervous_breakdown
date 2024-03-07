#!/bin/bash

## Author: Alicia Garcia-Roldan
## Date: March 2024

# This script is designed to rename bins files. It is perfect for SqueezeMeta (v1.6.2post3) SQM_project/results/bins folder. 
# You will need a directory.txt (with the names of the directory) and a samples.txt (with the names of the samples or the names you want to put as a prefix in your new bins)


DIRECTORY=$1
LINENUMBER=1

for DIRECTORY in $(cat directory.txt)
do

    echo -e  "\n It's Concha, I'm coming in (Soy Concha, entro)...${DIRECTORY} \n"

    SAMPLE=$(sed -n "${LINENUMBER}p" samples.txt)

    cd ${DIRECTORY}

    echo -e "\n You go around acting mystical, but you're bad, Yerbas, you're a creep \n (Vas de mística, pero eres mala Yerbas, eres un bicho) \n Prefix bin: ${SAMPLE} \n "
	
    prename "s/(.*.fa)/${SAMPLE}_\1/" *.fa
    prename "s/(.*.fa)/${SAMPLE}_\1/" *.fa.tax

    echo -e "I'm leaving..."
    
    cd ..

    LINENUMBER=$((LINENUMBER + 1))

done
