# somiteOrganoids
This repository contains code in three sections:
  1- Mathematical Model
  2- Quantification of Phase Profile on Anteroposterior Axis from Kymographs
  3- Single Cell RNA Sequencing Data Analysis
1- To run the simulation, you just need to run the MATLAB code.
2- For quantification of the phase profile, a sample Kymograph is provided with the code. The code takes the Kymograph as input and gives the anteroposterior phase profile as output.
3- Single cell sequencing analysis consists of two parts. First part is the normalization and clustering of the dataset. Second part is the pseudotime ordering of the analyzed dataset. To run the clustering code, download the dataset from NCBI GEO: GSE220563 and copy the files to the directory of the code. The first code will output adata.h5ad. The second code, analysis, takes adata.h5ad as input. 
