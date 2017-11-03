
############
Written by David Gfeller

The script requires MixMHCp to be installed in your path (see https://github.com/GfellerLab/MixMHCp)

For any question, please contact david.gfeller@unil.ch

MHCpExt can be used freely by academic groups for non-commercial purposes (see license).
The product is provided free of charge, and, therefore, on an "as is"
basis, without warranty of any kind.

FOR-PROFIT USERS
If you plan to use MHCpExt in any for-profit
application, you are required to obtain a separate license.
To do so, please contact eauffarth@licr.org or lfoit@licr.org at the Ludwig Institute for  Cancer Research Ltd.

If you use MixMHCpred1.1 in a publication, please cite:
Guillaume et al. The C-terminal extension landscape of naturally
presented HLA-I ligands, BioRxiv, doi:10.1101/213264 (2017).
Copyright (2017) David Gfeller

############


### Explanations ####

This scripts allows you to find cases of C- or N-terminal extensions in large datasets of peptide ligands, assuming that the canonical peptides have length 9.
It is tehreforeespecially designed to analyse HLA peptidomics datasets

As input, you should give a file containing a list of all peptides (any length). Typically, this could be the results of an HLA peptidomics screen on a cell line and tissue sample.

The script will first run MixMHCp on the 9-mers, providing an unsupervised view of the different motifs.
It will then identify longer peptides (10-, 11-, 12-mers) that cannot be explained by the bulge model derived from the 9-mers, but fit a N- or C-terminal extension model, and assign statisical significance (Z-score) to these cases.

To know which allele corresponds to which cluster, and to optimize the number of motifs, it is best to have a look at the logos found by MixMHCp and compare them with known logos.
Alternatively, you can run MixMHCp on the 9-mers and see in KLD/best_ncl.txt file the "optimal" number of motifs (this is the default value if the number of cluster is not specificied in MHCpExt).


### Installation ###

In MHCpExt, manually put the path to the MHCpExt/lib folder at "YOUR PATH TO MHCpExt/lib FOLDER".

Make sure awk and perl are installed on your computer and MixMHCp is installed and in your $PATH (see https://github.com/GfellerLab/MixMHCp).

The script should work on Mac and Linux. Depending on the size of the dataset, it can take several minutes to run.

### Running ###

Usage: MHCpExt -i INPUT_FILE -o OUTPUT_DIR

Mandatory parameters:

  -i, --input             Absolute or relative path to the input file (list of peptides)
  -o, --output            Name of the output directopry

Optional parameters: 

  -n, --ncl               Number of clusters to be used in MixMHCp - if not provided, the algorithm chooses it between 1 and 6
  -nm, --name             Name of the different clusters (separated by commas, e.g., A0101,A0201,B0702,B0801). Only valid if the number of clusters(-n) is given in input.


### Testing ###


To test your installation, from the MHCpExt/ directory, run:

MHCpExt -i test/A0301.txt -o test/out -n 1

This should take between one and two minutes on a standard computer.
Apart from path names, the output should be the same as in test/out_compare.
This file corresponds to peptides identified in HLA peptidomics study of a mono-allelic HLA-A03:01 positive cell line [Abelin et al Immunity 2017].
