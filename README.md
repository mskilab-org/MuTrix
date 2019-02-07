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

```{R}

muTrix                 package:muTrix                  R Documentation

muTrix

Description:

     Plots a matrix of mutations in genes x patients with pretty colors
     denoting mutation type.  Also works with seg data frames and
     dranger output to plot CNA's and rearrangements, respectively.

     Also works with lists of gene lists and lists of patients tracks
     (either named vectors, data frames, or matrices) to create panels
     of mutaitons.

     Assumes standard "annotated" maf file format (output of oncotator)
     with fields $Variant_Classification, $Hugo_Symbol

     Additional input is seg data frame with fields $ID, $chr, $start,
     $end corresponding to copy number segments for patient set.  If
     this is provided, a second panel will be drawn that contains
     amplified and deleted segme snts (ie with copy number greater /
     below cn.thresh / - cn.thresh) that are labeled broad vs focal
     based on focality.thresh.

Usage:

     muTrix(genes, patients = NULL, pat.tracks = NULL, bar.top = NULL,
       maf = NULL, mat = NULL, mut.categories = NULL,
       mat.colormap = topo.colors(20), dranger = NULL,
       pat.labels = F, srt.pat.labels = -90, cex.pat.labels = 0.9,
       cex.gene.labels = 1.1, include.blank = TRUE, keep.gene.ord = F,
       keep.pat.ord = F, mark.genes = NULL, plot.protein.change = FALSE,
       cex.protein.change = 1, srt.protein.change = 45, noncoding = TRUE,
       show.plot = TRUE, dump.protein.change = FALSE, cex.mut = 0.5,
       cex.cn = 0.95, cex.legend = 0.7, cex.tick = 1, cex.ylab = 1,
       srt.ylab = 90, adj.ylab = c(0.5, 0.5), 
       four.mut.categories = F, cex.pat.track = 1, pc.legend = 30,
       pc.sidebar = 10, pc.ylab = 30, freq.sidebar = T,
       top.marg.lim = NULL, bottom.marg.lim = NULL, single.scale = F,
       pad.track = 0.5, nchar.wrap = 20, pos.ylab = 0.3,
       gsub.ylab = c("[\\.\\_]", " "), lwd.border = 0.2,
       col.border = TRUE, col.canvas = "gray93",
       col.canvas.border = "gray50",
       col.pat.tracks = lapply(rownames(brewer.pal.info[brewer.pal.info$category
       == "qual", ]), function(x) brewer.pal(brewer.pal.info[x, "maxcolors"],
       x)), col.true = "darkgray", string.true = "True", ncol.legend = 2)
     
Arguments:

   genes: can be either (1) vector of genes (2) list of gene vectors,
          in which case each vector will be drawn as a separate panel
          with a label corresponding to the list name)

patients: optional character vector of patient ids specifying subset of
          patients to show

pat.tracks: optional named list of named vectors named by patient id
          specifying additional patient tracks or rownaems

     maf: "MAF" format data.frame with (at least) the columns
          $Tumor_Sample_Barcode, $Protein_Change, and
          $Variant_Classification, where Tumor_Sample_Barcode is the
          patient id

     mat: optional matrix of numeric data to plot, with columns named
          by patient ids

mat.colormap: vector mapping matrix values to colors

pat.labels: logical flag whether to show pat labels

srt.pat.labels: angle to rotate patient labels

plot.protein.change: flag whether to plot protein change

cex.protein.change: size of protein change label

srt.protein.change: rotation of protein change label

cex.legend: size of legend

mut.categories: optional named vector mapping variant classes (ie
          unique values in $Variant_Classification column of maf)to
          colors, used to specify alternate coloring schemes.

col.canvas: color of canvas

col.pat.tracks: list (same length as pat.tracks)
```
