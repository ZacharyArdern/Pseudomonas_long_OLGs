# Overview of data analysis scripts

# scripts are not complete programs but should be run by copying relevant sections into the Unix commandline 
# questions to: Zachary Ardern - zachary.ardern[]tum.de // z.ardern[]gmail.com - I am glad to help!

### Constructing and visualising gene trees: 

# find homologs of each mother gene with blastp and ncbi "IPG", align, create ML tree with IQtree
tree_and_alignment.bash 

# ete3 python script to visualise the alignment matched to the mother gene tree
visualise_tree_aln.py


### Stop codon analysis - "Frameshift" (using program from Schlub et al., 2018)

# script to run the R script below on each mother gene of interest, and to visualise results
Frameshift_analyses.bash

# modified R script; the script from Timothy Schlub is edited to incl. visualisation of simulated/permuted ORFs
Frameshift_20000_revcom0.r 


### Synonymous constraint analysis - "FRESCo" (using program from Sealfon et al., 2014)

# script to run FRESCo program on each mother gene of interest
FRESCo.bash


### OLG-appropriate dN/dS analysis - "OLGenie" (using program from Nelson et al., 2020)

# scripts to run FRESCo program on each mother gene of interest
OLGenie.bash

### Evolutionary simulation of 'mother' gene, observing alternate frame stops 
# (adapting a method from Cassan et al., 2016 - reimplemented in Pyvolve)

Simulation_analysis.bash

pyvolve_simulation_default-freqs.py 

### Analysis of protein evolution rates

protdist.bash