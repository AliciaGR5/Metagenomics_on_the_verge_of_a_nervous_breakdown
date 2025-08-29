# Metagenomics_on_the_verge_of_a_nervous_breakdown

<div align="center">
  <img width="1599" height="157" alt="image" src="https://github.com/user-attachments/assets/67af1113-cbdf-4d49-90ce-8104c3e61677" />
</div>

#### By Alicia Garcia-Roldan
:tropical_drink: Drink this _gazpacho_, and you'll understand metagenomics in a second :tropical_drink:

===========================================================================

**IMPORTANT INFORMATION:** These scripts are for:
  + SqueezeMeta v.1.6.2post3 
  + Bash v.4.4.20
  + R v.3.6.3 (2020-02-29) and v.4.1.2 (2021-11-01)

===========================================================================

**'Altruistic' cooperation among the prokaryotic community of Atlantic salterns assessed by metagenomics.**   
García-Roldán A, de la Haba RR, Sánchez-Porro C, Ventosa A. (2024). 
Microbiol Res. 288:127869. [doi: 10.1016/j.micres.2024.127869](https://www.sciencedirect.com/science/article/pii/S0944501324002702?via%3Dihub)<br><br>

:tropical_drink: You can see more scripts in the first part of this adventure!! Go and visit [**The_metagenomics_dispatch**](https://github.com/AliciaGR5/The_Metagenomics_dispatch) :tropical_drink:

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

Now you have your lovely bins...what can you do?? Maybe you should look for _important genes_...why? Maybe we find an important bin that does something cool, or maybe we find a relationship in the community...

_Remember: Nature is always singing "With a little help from my friends"_

Let's look for:

### [04. Biotin_synthesis.R](https://github.com/AliciaGR5/Metagenomics_on_the_verge_of_a_nervous_breakdown/blob/main/04.%20Biotin_synthesis.R)

Biotin synthesis genes

### [06. Biotin_extragenes.R](https://github.com/AliciaGR5/Metagenomics_on_the_verge_of_a_nervous_breakdown/blob/main/05.%20Biotin_extragenes.R)

Biotin synthesis genes and more...

### [07. Biotin_transporters.R](https://github.com/AliciaGR5/Metagenomics_on_the_verge_of_a_nervous_breakdown/blob/main/06.%20Biotin_transporters.R)

Biotin transporters genes

### [08. Bcarotene_synthesis.R](https://github.com/AliciaGR5/Metagenomics_on_the_verge_of_a_nervous_breakdown/blob/main/07.%20Bcarotene_synthesis.R)

β-carotene synthesis genes

### [09. Bactrub_synth.R](https://github.com/AliciaGR5/Metagenomics_on_the_verge_of_a_nervous_breakdown/blob/main/08.%20Bactrub_synth.R)

Bacterioruberin synthesis genes

### [10. Bcarot_rub_synth.R](https://github.com/AliciaGR5/Metagenomics_on_the_verge_of_a_nervous_breakdown/blob/main/09.%20Bcarot_rub_synth.R)

β-carotene and bacterioruberin synthesis genes

