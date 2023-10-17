#!/bin/bash

## Author: Alicia Garcia-Roldan
## Date: October 2023


FILE=$1

echo -e "I'm Concha, I'm coming in...(Soy Concha, entro)"

cd ~/tmp_bins

for FILE in $(cat ~/data/BINS/bad_bins_real.txt)
do
        echo -e "Come on, don't screw with me (Vamos, no me jodas)"
        rm ${FILE}

done

cd ..

echo -e "I'm Concha, I'm coming out...(Soy Concha, salgo)"

ls tmp_bins > after_rmv.txt

wc -l after_rmv.txt
wc -l ~/data/BINS/cool_bins.txt



# Before that, we should do

ls > everybin.txt # We obtain a list of the bins in the temp folder (those that have and don't hace quality values)
cut -d"." --complement -f5-  everybin.txt > everybin_var.txt # Delete the .f extension (we separate by dots, so we don't want the fifth column)

grep "^22IC" allbins_var_SQM.checkM > allbins_names.txt # We obtain a list with the bins from table 17 of SqueezeMeta (it was a checkm file so we only need the names)
cut -d" " --complement -f2-  allbins_name.txt > allbins_name_var.txt # Only, only the name, nothing more

# Parseamos

grep -f allbins_name_var.txt everybin_var.txt > cool_bins.txt # List with the coincidence
grep -vf allbins_name_var.txt everybin_var.txt > bad_bins.txt # List with the NO coincidence


sed 's/$/.fa/' bad_bins.txt > bad_bins_real.txt # We have to add again the .fa extension so that the script can find the bins that we want to remove


# Count the lines of the file of the bins that we DO want and the output of the script Rmv_bad_bins.sh...they should have the same number of lines

wc -l cool_bins.txt
wc -l allbins_name_var.txt