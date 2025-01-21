# Metagenomics_on_the_verge_of_a_nervous_breakdown
#### By Alicia Garcia-Roldan
:tropical_drink: Drink this gazpacho, and you'll understand metagenomics in a second :tropical_drink:

===========================================================================

**IMPORTANT INFORMATION:** These scripts are for:
  + SqueezeMeta v.1.6.2post3 
  + Bash v.4.4.20

===========================================================================

### [01. Change_name.sh](https://github.com/AliciaGR5/Metagenomics_on_the_verge_of_a_nervous_breakdown/blob/main/01.%20Change_name.sh)
This script is designed to rename the bins in a CheckM file. It is perfect for SqueezeMeta (v1.6.2post3) table 17.*.checkM located in intermediate folder.
Now you can concatenate all the files you want without having repeated names for different bins. 
You will need a files.txt (with the files you want to change) and samples.txt (with the names of the samples or the names you want to put in your new bins).

### [02. Rmv_bad_bins.sh](https://github.com/AliciaGR5/Metagenomics_on_the_verge_of_a_nervous_breakdown/blob/main/02.%20Rmv_bad_bins.sh)
This script is designed to play with lists of bins and write the perfect file for working with [mOTUlizer](https://github.com/moritzbuck/mOTUlizer).

### [03. Prefix_name.sh](https://github.com/AliciaGR5/Metagenomics_on_the_verge_of_a_nervous_breakdown/blob/main/03.%20Prefix_name.sh)
How do you play with bins from different metagenomes? Tired of making mistakes because the bins have the same name? Here's the solution!

This script is designed to change the name of all your bins so you can organise them by sample.

It is perfect for SqueezeMeta (v1.6.2post3) SQM_project/results/bins folder.

You will need a directory.txt (with the names of the directory) and a samples.txt (with the names of the samples or the names you want to put as a prefix in your new bins)

===========================================================================


