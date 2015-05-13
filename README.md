
#### Description

This repository contains supplementary materials for the manuscript titled "On the clustering of rare codons and its effect on translation" submitted to *In Silico Biology*.


#### Instructions

You must first download the compressed code (Code.zip) and data (Data.zip) files and extract them. 
The code has been tested using MATLAB version 7.6.0 (R2008a) on a Windows XP machine. 

- evalEcoli.m : this is the main script that detects clusters in Ecoli genes
- evalclust.m : this script evaluates various features of the detected clusters
- kscanstat.m, kscanstatmax.m, kscanstatmc.m : these functions perform the actual clustering procedure

You need not run the evalEcoli.m script, since all of its outputs have been provided as data files. 
To verify the results presented in the paper, simply run the script evalclust.m
