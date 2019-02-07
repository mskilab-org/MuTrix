<style type="text/css">
</style>

MuTrix
------

MuTrix (i.e. "Mutation-Matrix") allows you to render "oncoprint"-like
plots of mutation and clinical data from `data.frames` or `data.tables`
of `.maf` files and associated clinical data.

Note: Other packages such as
[ComplexHeatmap](https://jokergoo.github.io/ComplexHeatmap-reference/book/oncoprint.html)
are much more feature rich and likely more appropriate for modern cancer
genome analysis.


## Installation
------

1. Install devtools

```{r}
install.packages('devtools')
install.packages('testthat')
```
2. Install MuTrix

```{r}
devtools::install_github('mskilab/MuTrix)
```

## Usage
------

Example code to get you started:

```{R}

    library(data.table)
    library(muTrix)

    ## load small maf file
    maf = fread(system.file('extdata', "mutations.maf", package = "muTrix"))
    head(maf)
    ##    Tumor_Sample_Barcode Hugo_Symbol Variant_Classification Protein_Change
    ## 1:    LUAD-S01345-Tumor        BRD3      Missense_Mutation        p.Q531*
    ## 2:     LUAD-FH5PJ-Tumor         CBL      Missense_Mutation        p.G415S
    ## 3:    LUAD-S01345-Tumor         CBL      Missense_Mutation        p.Q175*
    ## 4:    LUAD-D02326-Tumor       FBXW7      Missense_Mutation        p.W365*
    ## 5:    LUAD-E01317-Tumor       FBXW7      Missense_Mutation        p.V464E
    ## 6:    LUAD-S01345-Tumor       FGFR3      Missense_Mutation        p.R750C
    ##    Variant_Type          id
    ## 1:          SNP LUAD-S01345
    ## 2:          SNP  LUAD-FH5PJ
    ## 3:          SNP LUAD-S01345
    ## 4:          SNP LUAD-D02326
    ## 5:          SNP LUAD-E01317
    ## 6:          SNP LUAD-S01345

    ## load clinical data, where "$id" matched the $Tumor_Sample_Barcode field in maf
    clinical = fread(system.file('extdata', "clinical.txt", package = "muTrix"))
    head(clinical)
    ##                   id smoking_status stage gender
    ## 1: LUAD-S01345-Tumor      Ex Smoker    IB   Male
    ## 2:  LUAD-FH5PJ-Tumor      Ex Smoker    IB   Male
    ## 3: LUAD-D02326-Tumor   Never Smoker   IIA Female
    ## 4: LUAD-E01317-Tumor      Ex Smoker       Female
    ## 5: LUAD-S01302-Tumor Current Smoker    IB   Male
    ## 6: LUAD-S00488-Tumor      Ex Smoker  IIIA   Male

    ## pick genes
    genelists = list(
              Oncogenes = c('EGFR', 'KRAS', 'PIK3CA'),
              TSG = c('TP53', 'STK11', 'NF1'),
              Other = c('KEAP1', 'RBM10', 'U2AF1', 'SMARCA4'))
    genelists
    ## $Oncogenes
    ## [1] "EGFR"   "KRAS"   "PIK3CA"
    ## 
    ## $TSG
    ## [1] "TP53"  "STK11" "NF1"  
    ## 
    ## $Other
    ## [1] "KEAP1"   "RBM10"   "U2AF1"   "SMARCA4"

    ## example muTrix plot over above maf, genelists, and patient metadata
    mat = muTrix(maf = maf,
           genes = genelists,
           pat.tracks = list(
                 gender = clinical[, structure(gender, names = id)], ## vectors named by patient id
                 smoking = clinical[, structure(smoking_status, names = id)]))
```

![](https://raw.githubusercontent.com/mskilab/MuTrix/master/tutorial_files/figure-markdown_strict/unnamed-chunk-5-1.png)

```{R}
## for more documentation
?muTrix
```
