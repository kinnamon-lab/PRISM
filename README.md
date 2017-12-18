# PRISM: Polygenic RIsk Score Modeler
Copyright 2014-2017 The Ohio State University Wexner Medical Center\
**Authors:** Daniel D. Kinnamon, PhD; Carl A. Starkey, PhD\
Licensed under the Apache License, Version 2.0

## INTRODUCTION

PRISM is a set of Java classes that can be used to implement polygenic risk score models given a set of SNPs and an annual incidence file for a disease in the same population. In addition to providing classes useful for building such risk models, it provides an executable `RiskModelBuilder` class for fitting new risk model objects and another executable `RiskPredictor` class for generating risk predictions for individuals using these risk model objects. Additional details are provided below.

## LICENSE

See [here](LICENSE.md). This software is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied.

## CITATION

If you use the PRISM software in your research, please include the following citation:

    Kinnamon DD, Starkey CA. PRISM: Polygenic RIsk Score Modeler [Source code on Internet]. Version X. Columbus (OH): The Ohio State University Wexner Medical Center; DD MMM YYYY version release date [cited DD MMM YYYY download date]. Available from: https://github.com/kinnamon-lab/PRISM.

## INSTALLATION
Coming soon.

## RUNNING `RiskModelBuilder`

### USAGE

A detailed usage message for the `RiskModelBuilder` class can be obtained by running the following on the command line (Linux shells or MS-DOS):

    java -cp /path/to/prism.jar edu.osumc.prism.RiskModelBuilder -h

The RiskModelBuilder needs two types of input files to build a risk model: a SNP information file and an annual incidence file. These files should be named `[model ID]_SNPs.dat` and `[model_ID]_annInc.dat`, where `[model ID]` is the name that you would like the fitted model and risk model object to have. A successful run of `RiskModelBuilder` will produce one or more serialized risk model object files, which are binary files that have the extension `*.rmo`. These contain all information necessary to generate relative and absolute risk predictions for one or more individuals using the `RiskPredictor` class.

### SNP INFORMATION FILE

The `[model ID]_SNPs.dat` file should contain the following columns in the specified order:

Column Name | Data Type | Description
----------- | --------- | -----------
modelID | String | The model identifier used in the file name
rsID | String | The dbSNP RefSNP identifier. Must match regular expression `^rs[0-9]+$`
sourcePub | String | The publication or other source supplying log HR and allele frequency data for the SNP
allele1 | String | Allele 1 in the source. See below for additional information
allele2 | String | Allele 2 in the source. See below for additional information
orientRs | String | Orientation of the alleles in the source relative to refSNP alleles. See below for additional information
allele2Freq | Double | *Allele 2* frequency from source publication or a public database
allele2LnHR | Double | Cox proportional hazards model natural log hazard ratio *for each additional allele 2* obtained from the source

The names of these columns should appear in a single header row followed by one row for each SNP for a particular model.

Alleles should be uppercase character strings comprising one or more A, C, G, or T characters, a single - for a deletion, or a single 0 for a missing allele. Only biallelic SNPs are supported in risk models. Orientation relative to RefSNP alleles is "Forward" if allele calls are on the same DNA strand as the RefSNP alleles. It is "Reverse" if they are on the opposite strand. This orientation should not be confused with orientation relative to the reference sequence, contig, or anything else, which all may differ. For example, take rs2420946 with RefSNP alleles C/T in forward orientation relative to the reference sequence (see dbSNP). If the source provides C/T as the possible alleles, the reported data would be in forward orientation relative to the RefSNP alleles. If the source provides G/A as the possible alleles, the reported data would be in reverse orientation. Now take rs1292011, which also has RefSNP alleles C/T but in reverse orientation relative to the reference sequence (see dbSNP). If the source provided C/T as the possible alleles, *the reported data would still be in forward orientation relative to the RefSNP alleles*. If the source provided G/A as the possible alleles, the reported data would still be reverse orientation.

### ANNUAL INCIDENCE FILE

The `[model ID]_annInc.dat` file should contain the following columns in the specified order:

Column Name | Data Type | Description
----------- | --------- | -----------
modelID | String | The model identifier used in the file name
ageYrs | Integer | Age in years
annInc | Double | The annual incidence between ages ageYrs - 1 and ageYrs. This must be 0 for ageYrs = 0

The names of these columns should appear in a single header row followed by one row for each consecutive year between 0 and the last age for which predictors are desired, inclusive

## RUNNING `RiskPredictor`

### USAGE

A detailed usage message for the `RiskPredictor` class can be obtained by running the following on the command line (Linux shells or MS-DOS):

    java -cp /path/to/prism.jar edu.osumc.prism.RiskPredictor -h

The RiskPredictor needs three input files to produce risk predictions: a risk model object file `[model ID].rmo` generated by the `RiskModelBuilder` class, a PED file `[name].ped` containing genotypes for the SNPs in this risk model, and a MAP file `[name].map` describing the SNPs in the PED file. Running `RiskPredictor` on these files will generate an output file `[name]-[model ID].prd` containing risk predictions for the individuals in the PED file.

### PED FILE
The PED file should be a tab-delimited file containing no header row, 1 person per row, and the following columns in the specified order:

Column | Data Type | Description
----------- | --------- | -----------
Individual ID | String | Individual IDs are alphanumeric and should uniquely identify a person
SNP 1, Allele 1 | String | Allele 1 for individual at first SNP in MAP file
SNP 1, Allele 2 | String | Allele 2 for individual at first SNP in MAP file
... | ... | ...
SNP k, Allele 1 | String | Allele 1 for individual at kth SNP in MAP file
SNP k, Allele 2 | String | Allele 2 for individual at kth SNP in MAP file

The same restrictions on alleles in the SNP information file apply here. All SNPs must have two alleles specified. Either both alleles should be missing (i.e. 0) or both should be supplied. IUPAC codes for ambiguous nucleotides may not be used; these genotypes must be set to missing. An example PED file containing two individuals typed for 3 SNPs is included below: 

    OSU1234	A	A	G	G	A	C
    MOF1234	A	A	-	G	0	0

### MAP FILE
The MAP file should be a tab-delimited file containing no header row, one SNP per row (in the same order as the genotype columns of the PED file are read from left to right), and the following columns:

Column | Data Type | Description
----------- | --------- | -----------
rs # | String | dbSNP RefSNP ID
orientRs | String | Orientation of assay alleles relative to RefSNP alleles: "Forward" or "Reverse"

An example MAP file corresponding to the above PED file is included below:

    rs123456	Forward
    rs234567	Reverse

### OUTPUT FILE

The output file for a particular PED/MAP pair and `*.rmo` object is a tab-delimited file with a header row and one row per individual and the following columns:

Column Name | Data Type | Description
----------- | --------- | -----------
IndivID | String | ID from first column of PED file
Model | String | Name of the risk model
rs###### (multiple) | String | Genotype at this SNP
PI | Double | Polygenic risk score (Cox linear predictor from risk model)
PIPctl | Double | Estimated percentile of this polygenic risk score in the population assuming Hardy-Weinberg equilibrium at and linkage equilibrium between all model SNPs
PredCumRiskAge## (multiple) | Double | Predicted cumulative risk (1 - survivor function estimate) for this individual at age ##

## METHODOLOGICAL DETAILS
Coming soon.

## ACKNOWLEDGEMENTS

Development of this software was supported in part by Award Number Grant UL1TR001070 from the National Center For Advancing Translational Sciences. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Center For Advancing Translational Sciences or the National Institutes of Health.