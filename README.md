# EMMAX Pipeline
This repo includes a number of scripts used to perform a mixed model regression GWAS, with alternative model masking, as published in Ravenhall *et al*. 2018 https://doi.org/10.1371/journal.pgen.1007172. It is primarily meant as a personal archive of those scripts for later use.

## Core Pipeline
### EMMAX_Pipeline.sh
- Core pipeline for performing EMMAX mixed model regression with allele masking.
- Input:
    - `CHROMOSOME` - Chromosome, first command-line argument
    - `SUBTYPE` - Subtype, second command-line argument
    - `EMMAXdir` - location of emmax
    - `TPED` - in 12 format
    - `PHENO` - phenotype file: <ID> <ID> <CASE/CONTROL>, no header
    - `COVARS` - covariates file: <ID> <ID> <CASE/CONTROL> <COVARIATES>, no header
    - `OUTDIR` - output directory
    - `OUTPREFIX` - output prefix
    - `KINF` - kinship matrix
    - `SCRIPTDIR` - directory with supporting scripts

### mergePlotEMMAX.r
- Manhattan plot based on EMMAX output
- Supporting script for EMMAX_Pipeline.sh

## Associated Scripts
### plotSNPIntensity.r
- Plot SNP intensity for a given SNP, useful for validating classification for key candidates.
- Input:
    - `out.XY.FORMAT`: Rows as SNPs, Columns as 'CHR:Chromosome, POS:Basepair, ID_0:Sample0, ..., ID_N:SampleN', calls as 'intensityX,intensityY'.
    - `out.GT.FORMAT`: Rows as SNPs, Columns as 'CHR:Chromosome, POS:Basepair, ID_0:Sample0, ..., ID_N:SampleN', calls as ./., 0/1, 1/1 etc.
    - `PLACEHOLDER.sample`: Case/Control sample file, requires 'scanID' (IDs) and 'caseorcontrol' ('CASE'/'CONTROL') columns

### getLambda.r
- Scrapbook of functions for calculating genomic inflation factor (lambda).
- Requires a .ps file as input.
- Functionality is present within mergePlotEMMAX.r

### ProduceTopSNPs.r
- Create a refined subset of the most significant SNPs from a combined EMMAX pipeline output.

### runImpute2.sh
- Run impute2 in parallel
- Input:
    - `impute2dir` - location of impute2
    - `refDir` - reference file directory, containing genetic map, .hap and .legend files 
    - `refName` - prefix of reference files
    - `mapNamePre` - Pre-chromosome component of genetic map file name
    - `mapNamePost` - Post-chromosome component of genetic map file name
    - `toImpFile` - .gen file for inputation
    - `outprefix` - output prefix
    - `psLimit` - number of sub-processes to spawn
