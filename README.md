===============================================================================
                      PRISM - Polygenic RIsk Score Models
                               * Version 0.2.5 *
     -- Copyright 2014 The Ohio State University Wexner Medical Center --
                Licensed under the Apache License, Version 2.0
===============================================================================
NOTE: This software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied.
-------------------------------------------------------------------------------

usage: java -jar prism-0.2.5.jar -h | -m | -p <fileroot>
 -h,--help                 Print this usage information.
 -m,--models               Print available risk model information to
                           console.
 -p,--predict <fileroot>   Generate risk predictions for individual
                           genotype data supplied in the PED/MAP file pair
                           given by <fileroot>.ped/<fileroot>.map.
                           Predictions will be output to
                           <fileroot>-<modelID>.prd files, one for each
                           risk model, and logging information will be
                           printed to the console.

--------------------------------------------------------------------------------

PED file
----------
The PED file should be a tab-delimited file containing no header row, 1 person 
per row, and the following columns:
    - Individual ID: Individual IDs are alphanumeric and should uniquely 
      identify a person.
    - Genotypes (column 2 onwards): Should also be tab-delimited with two 
      columns per SNP (one for each allele). The alleles should be character 
      strings (e.g, A, C, G, T, -) or 0 (the missing genotype character). All 
      markers should be biallelic, and all SNPs must have two alleles specified.
      Either both alleles should be missing (i.e. 0) or neither.

As an example, here is what a file containing two individuals typed for 3 SNPs 
would look like:
OSU1234 A A G G A C
MOF1234 A A - G 0 0

IMPORTANT: We cannot use IUPAC codes for ambiguous nucleotides. If the assay 
doesn't produce them, then there's no problem, but those genotypes will have to 
be set to missing otherwise.

MAP file
----------
The MAP file should be a tab-delimited file containing no header row, one marker
per row (in the same order as the genotype columns of the PED file are read 
from left to right), and the following columns:
    - rs#
    - Orientation relative to refSNP in assay: "Forward" or "Reverse" (see 
      below)

As an example, here is what the MAP file corresponding to the above PED file 
might look like:
rs123456 Forward
rs234567 Reverse

IMPORTANT: Orientation relative to RefSNP alleles is "Forward" if your assay's 
allele calls are on the same DNA strand as the RefSNP alleles. It is "Reverse" 
if they are on the opposite strand. This orientation should not be confused with
orientation relative to the reference sequence, contig, or anything else, which 
all may differ. For example, take rs2420946 with RefSNP alleles C/T in forward 
orientation relative to the reference sequence (see dbSNP). If your assay called
C/T as the possible alleles, it would be in forward orientation relative to the 
RefSNP alleles. If it called G/A as the possible alleles, it would be in reverse
orientation. Now take rs1292011, which also has RefSNP alleles C/T but in 
reverse orientation relative to the reference sequence (see dbSNP). If your 
assay called C/T as the possible alleles, **it would still be in forward 
orientation relative to the RefSNP alleles**. If it called G/A as the possible 
alleles, it would still be reverse orientation.
