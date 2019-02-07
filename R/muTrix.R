#' @name muTrix
#' @title muTrix
#'
#' @description
#' Plots a matrix of mutations in genes x patients with pretty colors denoting mutation type.  Also works with seg data frames and dranger output to plot CNA's and rearrangements, respectively. 
#'
#' Also works with lists of gene lists and lists of patients tracks (either named vectors, data frames, or matrices) to create panels of mutaitons.
#'
#' Assumes standard "annotated" maf file format (output of oncotator) with fields $Variant_Classification, $Hugo_Symbol
#'
#' Additional input is seg data frame with fields $ID, $chr, $start, $end corresponding to copy number segments for patient set.  If this is provided, a second
#' panel will be drawn that contains amplified and deleted segme snts (ie with copy number greater / below cn.thresh / - cn.thresh) that are labeled
#' broad vs focal based on focality.thresh.  
#'
#' @param maf "MAF" format data.frame with (at least) the columns $Tumor_Sample_Barcode, $Protein_Change, and $Variant_Classification, where Tumor_Sample_Barcode is the patient id
#' @param genes can be either (1) vector of genes  (2) list of gene vectors, in which case each vector will be drawn as a separate panel with a label corresponding to the list name)
#' @param patients optional character vector of patient ids specifying subset of patients to show
#' @param mut.categories optional named vector mapping variant classes (ie unique values in $Variant_Classification column of maf)to colors, used to specify alternate coloring schemes.
#' @param pat.tracks optional named list of named vectors named by patient id specifying additional patient tracks or rownaems
#' @param mat optional matrix of numeric data to plot, with columns named by patient ids
#' @param mat.colormap vector mapping matrix values to colors
#' @param col.pat.tracks list (same length as pat.tracks)
#' @param pat.labels logical flag whether to show pat labels
#' @param col.canvas color of canvas
#' @param srt.pat.labels angle to rotate patient labels
#' @param plot.protein.change flag whether to plot protein change
#' @param cex.protein.change size of protein change label
#' @param srt.protein.change rotation of protein change label
#' @param cex.legend size of legend
#' @return matrrix
muTrix = function(genes, # can be either (1) vector of genes  (2) list of gene vectors
                         # (in which case each vector will be drawn as a separate panel with a label corresponding to the list name)
  patients = NULL,
  pat.tracks = NULL, # this is either
                     # (1) a vector of patient categories whose names are patient label
                     # (2) a matrix of features x patients or
                     # (3) a list of vectors or matrices (in which case there will be one panel drawn per list, labeled by the names
  bar.top = NULL, # bar plot data to put on top - this is a list with fields $data (matrix, df, or named vector of data x patients), names have to match ID's in maf and seg
  maf = NULL,
  seg = NULL,
  mat = NULL, ## override for seg [this arg ignored if seg is provided] with two possibilities:
              ## (1) named numeric matrix of genes x patients continuous values, which will be plotted in background by interpolating mat.colormap from lowest (first) to highest (last) value (color)
              ## (2) named character matrix of genes x patients continuous values, which will be plotted in background as per (named) mat.colormap vector                     
  mat.colormap = topo.colors(20), ## if mat is character names of mut.colormap should match values of mat, otherwise can be unnamed vector of colors that will be interpolated from lowest to highest
  seg.field = 'seg.mean', cn.thresh = 0.2, cn.thresh.del = -cn.thresh, cn.thresh.amp = cn.thresh, focality.thresh = 5e6, # seg processing options
  rg = if (!is.null(seg)) read_refGene(hg19 = T), # used for mapping segs to genes, can be preloaded for faster use (RECOMMENDED) .. default assumes hg19 segs
  pat.labels = F, srt.pat.labels = -90, cex.pat.labels = 0.9, cex.gene.labels = 1.1, include.blank = TRUE,
  mut.categories = NULL,  # named vector of colors or un-named vector of categories, determines how to interpret values in Variant Classification field
  keep.gene.ord = F, keep.pat.ord = F, # keep supplied order, instead of arranging genes and patients in staircase pattern
  mark.genes = NULL, # either (1) vector of genes which to mark (ie with a star symbol) or (2) vector of R points (i.e. pch) named by genes telling which R symbol to mark which gene with  
  plot.protein.change = FALSE, # will plot protein change on top of colors (won't look pretty unless very few patients)
  cex.protein.change = 1,
  srt.protein.change = 45,
  dranger = NULL, # UNDER CONSTRUCTION CURRENTLY!
  noncoding = TRUE,     
  show.plot = TRUE,  # if FALSE, then will only output the mutrix struct and not draw any plot 
  dump.protein.change = FALSE, # will dump protein change (e.g. p.G12V) instead of variant classification
  cex.mut = 0.5, # height of mutation square (should be smaller than cn and < 1)
  cex.cn = 0.95, # height of cn rectangle (should be larger than mut and < 1)
  cex.legend = 0.7, # size of legend
  cex.tick = 1, # tick sizes for various scales on plot 
  cex.ylab = 1, # ylabel size (ie labeling each submatrix)  
  srt.ylab = 90, # ylabel orientation
  adj.ylab = c(0.5, 0.5), # ylabel adj
  four.mut.categories = F, # if true, will only output silent, missense, truncating, in-frame indel
  cex.pat.track = 1, # size which patient tracks should be scaled 
  pc.legend = 30, # percentage of the horizontal plot real estate that the legend takes up (must be a number between 0 and 100)
  pc.sidebar = 10, # percentage of horizontal plot occupied by sidebar stacked bargraph (on right of main matrix)
  pc.ylab = 30, # pc of horizontal real estate left to the gene / patient track labels
  freq.sidebar = T, # if T then will plot marginal counts as % of cohort, otherwise will plot as raw numbers, 
  top.marg.lim = NULL, # side bar x limits corresponding to marginal counts of unique values of mutation tracks
  bottom.marg.lim = NULL, # side bar x limits corresponding to marginal counts of unique values of patient tracks
  single.scale = F, # will scale top and bottom marg limits to be the same (ie the max of the two)
  pad.track = 0.5, # number of rows of whitespace to insert between each track (e.g. genelist, patient track)
  nchar.wrap = 20, # number of characters to wrap around text (i.e. y labels)
  pos.ylab = 0.3, # ylabel position (relative to matrix width)
  gsub.ylab = c('[\\.\\_]', ' '), # gsub first item with second item in ylabels
  lwd.border = 0.2, # width of border around mutation square
  col.border = TRUE,
  col.canvas = 'gray93',
  col.canvas.border = 'gray50', 
  col.pat.tracks = lapply(rownames(brewer.pal.info[brewer.pal.info$category == "qual", ]), function(x) brewer.pal(brewer.pal.info[x, 'maxcolors'], x)), # (alternating) list of color palettes to use for (non boolean) patient tracks
  col.true = 'darkgray', # color to use for "true" entries in boolean patient tracks
  string.true = 'True',
  ncol.legend = 2)  # string to use in legend for "true" entries in boolean patient tracks , 
  {
    suppressPackageStartupMessages(require(RColorBrewer))
    suppressPackageStartupMessages(require(gplots))

    if (is.null(genes) | length(genes)==0)
      stop('Empty gene list')
    
    if (is.list(genes)) # if genes provide as list then we will plot each separately
      {
        gene.lists = genes;
        genes = unique(unlist(genes))
      }
    else # otherwise we have a single (trivial) list
      gene.lists = list(genes)
    
    genes = setdiff(genes, '');
    genes.og = genes;
    
    ## TEMPORARY
    if (!is.null(dranger))
      {
        warning('Rearrangement plotting currently under construction .. ignoring dranger input for now. ')
        dranger = NULL;
      }
    
    if (is.null(maf) & is.null(seg) & is.null(mat))
      stop('Either maf or seg or mat must be specified!');

    #########
    # BUILD variant matrices
    # for muts, rearrangements, etc
    #########    
    if (!is.null(maf))
     {

         if (!is.data.frame(maf))
             maf = as.data.frame(maf)
                                        # get patient name (different naming conventions in different maf files)
        if (is.null(maf$patient_name))
          if (is.null(maf$name))
            if (is.null(maf$patient))
                maf$patient_name = maf$Tumor_Sample_Barcode
#              maf$patient_name = gsub('\\-Tumor', '', maf$Tumor_Sample_Barcode)
            else
              maf$patient_name = maf$patient
          else
            maf$patient_name = maf$name;

        maf = maf[!is.na(maf$patient_name), ]
        maf.subset = maf[maf$Hugo_Symbol %in% genes, , drop = FALSE]

                                       # record mutation state numeric codes (will replace back in later)        
        mut.map = list()

        if (!is.null(mut.categories))
          {
            if (is.null(names(mut.categories)))
              mut.categories = structure(names = mut.categories, brewer.master(length(mut.categories)))

            mut.map[as.character(c(0:length(mut.categories)))] = c('', names(mut.categories));            
            MIN.NS = 2;

            mut.var.colors = mut.categories;
          }        
        else if (four.mut.categories)
          {
            mut.map[as.character(c(0:4))] = c("", "Silent", "In Frame Indel", "Missense", "Truncating");
            tmp <- lapply(list(c("_Mutation", ""), c("Silent", "Silent"), c("Nonsense", "Truncating"),
                               c("Stop_Codon", "Truncating"),
                               c("Start_Codon", "Truncating"),
                               c("Splice_Site.*", "Truncating"), 
                               c("Nonstop|De_novo_Start_.*|Start_.*|.*Flank", ""),
                               c("[iI]n_[fF]rame.*", "In Frame Indel"),
                               c("Frame_Shift.*", "Truncating"), c('.*(UTR)|(Intron).*', '')),
                          function(x) maf.subset$Variant_Classification <<- sub(x[1], x[2], as.character(maf.subset$Variant_Classification)))
            MIN.NS = 2;

            mut.var.colors = c(lighten(brewer.pal(9, "Set1")[c(3,2,1,5)], 0))
            names(mut.var.colors) = c("Silent", "Missense", "Truncating", "In Frame Indel");
            mut.var.colors['In Frame Indel'] = 'cyan'
          }
        else
          {
            mut.map[as.character(c(0:8))] = c("", "Silent.",  "Noncoding",  "In Frame Indel",
                     "Other non syn.", "Missense", "Splice Site", "Frame Shift", "Nonsense");
            tmp <- lapply(list(c("_Mutation", ""), c("Silent", "Syn."), c("Splice_Site.*", "Splice Site"), 
                               c("Nonstop|De_novo_Start_.*|Start_.*|Stop_.*", "Other non syn."), c("[iI]n_[fF]rame.*", "In Frame Indel"),
                           c("Frame_Shift.*", "Frame Shift"), c('.*((UTR)|(Intron)|(IGR)).*', 'Noncoding')),
                      function(x) maf.subset$Variant_Classification <<- sub(x[1], x[2], as.character(maf.subset$Variant_Classification)))
            MIN.NS = 3;

            mut.var.colors = c(lighten(brewer.pal(9, "Set1")[c(3,2,4,1,5,6,7)], 0))
            names(mut.var.colors) = c("Silent", "Missense", "Splice Site", "Nonsense", "Frame Shift", "In Frame Indel", "Other non syn.");
            mut.var.colors['In Frame Indel'] = 'cyan'
            mut.var.colors['Noncoding'] = 'gray35'
          }
        
        maf.subset$Variant_Classification_Num <- as.numeric(factor(maf.subset$Variant_Classification, 
                                                            levels=as.character(mut.map[as.character(1:(length(mut.map)-1))])))

                 
        maf.subset.agg <- aggregate(formula = Variant_Classification_Num ~ Hugo_Symbol + patient_name, FUN = max, na.rm=TRUE, data=maf.subset)
        genematrix.mut <- xtabs(formula = Variant_Classification_Num ~ Hugo_Symbol + patient_name, data=maf.subset.agg, drop.unused.levels = FALSE)

        if (is.null(patients))
          {
            patients = unique(maf$patient_name);

            if (!is.null(seg))
              patients = union(patients, seg$ID)
            else if (!is.null(mat))
              patients = union(patients, colnames(mat))

            if (!is.null(dranger))
              patients = union(patients, dranger$individual)
          }
        else
          patients = setdiff(patients, NA)
        
        others = setdiff(patients, maf.subset.agg$patient_name);
        genematrix.mut = cbind(genematrix.mut, matrix(0, nrow = nrow(genematrix.mut), ncol = length(others), dimnames = list(rownames(genematrix.mut), others)))
        genematrix.mut = genematrix.mut[, patients, drop = FALSE];
        
                                        # fill in "blank" genes
        others = setdiff(genes, rownames(genematrix.mut));
        genematrix.mut = rbind(genematrix.mut, matrix(0, ncol = ncol(genematrix.mut), nrow = length(others), dimnames = list(others, colnames(genematrix.mut))))
      }
    
    if (!is.null(seg))
      {
        REQUIRED.SEG.FIELDS = c('ID', 'chr', 'pos1', 'pos2');
                        
        seg = standardize_segs(seg)
        seg$length = seg$pos2-seg$pos1;

        for (required.field in c(REQUIRED.SEG.FIELDS))
          if (!(required.field %in% names(seg)))
            stop(sprintf('Error: Field "%s" not found in the columns of seg data frame.', required.field))

          if (!(seg.field %in% names(seg)))
            stop(sprintf('Error: Segment value field "%s" not found in the columns of seg data frame.\n\tPlease add this column to the seg data frame or specify an alternate seg.field. ', seg.field))

        if (is.null(patients))
          patients = union(union(maf$patient_name, seg$ID), dranger$individual)

        # include only segs matching these genes
        ix = match(genes, rg$gene_sym);
        genes.df = data.frame(gene = genes, stringsAsFactors = F);  
        genes.df[, c('chr', 'pos1', 'pos2')] = rg[ix, c('chr', 's1', 'e1')]
        
        seg2 = seg[seg.on.seg(seg, genes.df),];

        genematrix.broad.amp = cn_matrix(genes, seg[seg$length >= focality.thresh, ], agg = "max", patients = patients, field = seg.field,  rg = rg)>cn.thresh.amp
        genematrix.broad.del = cn_matrix(genes, seg[seg$length >= focality.thresh, ], agg = "min", patients = patients, field = seg.field, rg = rg)<cn.thresh.del
        genematrix.focal.amp = cn_matrix(genes, seg[seg$length < focality.thresh, ], agg = "max", patients = patients, field = seg.field, rg = rg)>cn.thresh.amp
        genematrix.focal.del = cn_matrix(genes, seg[seg$length < focality.thresh, ], agg = "min", patients = patients, field = seg.field, rg = rg)<cn.thresh.del
        genematrix.focal.na = is.na(cn_matrix(genes, seg, agg = "min", patients = patients, field = seg.field, rg = rg, use.na = TRUE))

        ## code copy number events
        genematrix.cn = genematrix.broad.amp;
        genematrix.cn[genematrix.broad.del] = 2
        genematrix.cn[genematrix.focal.amp] = 3
        genematrix.cn[genematrix.focal.del] = 4
        genematrix.cn[((genematrix.focal.amp + genematrix.broad.amp)*(genematrix.broad.del + genematrix.focal.del))>0] = 5
        genematrix.cn[,colSums(!genematrix.focal.na)==0] = NA;

        ## record copy number codes as list (will replace back in later)       
        cn.map = list();
        cn.map[as.character(c(0:5))] = c("", "Broad Amp", "Broad Del", "Focal Amp", "Focal Del", "Amp + Del");        
#        cn.map = list(c('0', ''), c('1', 'Broad Amp'), c('2', 'Broad Del'), c('3', 'Focal Amp'), c('4', 'Focal Del'), c('5', 'Amp+Del'))         
      }
    else if (!is.null(mat))
      {
        ## pad mat with any missing rows / cols
        if (length(x <- setdiff(genes, rownames(mat)))>0)
          mat = rbind(mat, matrix(NA, nrow = length(x), ncol = ncol(mat), dimnames = list(x, colnames(mat))))

        if (length(x <- setdiff(patients, colnames(mat)))>0)
          mat = cbind(mat, matrix(NA, ncol = length(x),  nrow = nrow(mat), dimnames = list(rownames(mat), x)))

        # enforce order
        mat = mat[genes, patients]
      }

    if (!is.null(dranger))
      {
        if (is.null(patients))
          patients = union(union(maf$patient_name, seg$ID), dranger$individual)
        
        dranger$Variant_Classification = gsub('([^\\(]+)[\\(]?.*$', '\\1', dranger$fusion)
        
                                        # record mutation state numeric codes (will replace back in later)
        tmp <- lapply(list(c("^-$", "Intergenic"), c("Deletion of .*", "In-frame del"), c("Duplication of .*", "In-frame dup"),
                           c('Protein fusion: in frame.*', 'In-frame fusion'), c('Protein fusion: out of frame.*', 'Out-of-frame fusion'),
                           c('Duplication within transcript: mid-exon.*', 'Mid-exon duplication'),
                           c('Deletion within transcript: mid-exon.*', 'Mid-exon deletion'),
                           c('Deletion within intron', 'Intronic del / dup'),
                           c('Deletion within transcript', 'In-frame del'),
                           c('Duplication within transcript', 'In-frame dup'),
                           c('Duplication within intron', 'Intronic del / dup'),
                           c('Protein fusion: mid-exon.*', 'In-frame fusion')),
                      function(x) dranger$Variant_Classification <<- sub(x[1], x[2], as.character(dranger$Variant_Classification)))

        ra.uvariants = c('Intergenic', 'Intronic del / dup', 'Mid-exon deletion', 'Mid-exon duplication','In-frame del', 'In-frame dup', 'Transcript fusion ', 'Antisense fusion', 'Out-of-frame fusion', 'In-frame fusion')
        ra.uvariants = c(setdiff(dranger$Variant_Classification, ra.uvariants), ra.uvariants)
        
        dranger.subset <- subset(dranger, dranger$gene1 %in% genes | dranger$gene2 %in% genes)
        dranger.subset$Variant_Classification = factor(dranger.subset$Variant_Classification, levels = ra.uvariants)
        dranger.subset$Variant_Classification_Num = as.numeric(dranger.subset$Variant_Classification)
        dranger.subset$label <- paste(dranger.subset$gene1, dranger.subset$gene2, sep = "-");
        dranger.subset$label[dranger.subset$gene1 == dranger.subset$gene2] = dranger.subset$gene1[dranger.subset$gene1 == dranger.subset$gene2];
        dranger.subset$label[is.na(dranger.subset$label)] = '';
        ra.ulabels = unique(dranger.subset$label);        
        dranger.subset$label <- factor(dranger.subset$label, levels=ra.ulabels)
        
        dranger.subset <- aggregate(Variant_Classification_Num ~ label + individual, max, na.rm=TRUE, data=dranger.subset) 
        genematrix.ra <- xtabs(Variant_Classification_Num ~ label + individual, data=dranger.subset, drop.unused.levels = FALSE)

        others = setdiff(patients, dranger.subset$individual);
        genematrix.ra = cbind(genematrix.ra, matrix(0, nrow = nrow(genematrix.ra), ncol = length(others), dimnames = list(rownames(genematrix.ra), others)))
        genematrix.ra = genematrix.ra[, patients, drop = FALSE];        
      }

    if (length(patients)==0)
      stop('No patients to plot!')
    
    
    #########
    # ORDER GENES and PATIENTS
    # 
    #########

    if (noncoding)
        MIN.NS = 1
    
    if (!is.null(maf))
      genes = union(genes, rownames(genematrix.mut)[rowSums(genematrix.mut)>0])
    
    if (!is.null(seg))
      genes = union(genes, rownames(genematrix.cn)[rowSums(genematrix.cn, na.rm = TRUE)!=0])
    else if (!is.null(mat))
      genes = union(genes, rownames(mat)[rowSums(!is.na(mat), na.rm = TRUE)])
          
    if (keep.gene.ord)
      gene.ord = intersect(genes.og, genes)
    else
      if (!is.null(maf))
        gene.ord = rownames(genematrix.mut)[order(-rowSums(genematrix.mut>=MIN.NS))] ## order using non-silent coding variants (hacky coding for now)
      else if (!is.null(seg))
        gene.ord = rownames(genematrix.cn)[order(-rowSums(genematrix.cn!=0))]
      else if (!is.null(mat))
        gene.ord = rownames(mat)[order(-rowSums(mat!=0))]
      else if (!is.null(dranger))
        gene.ord = rownames(genematrix.cn)[order(-rowSums(genematrix.ra!=0))]
    
    if (!include.blank)
      {
        tmp = NULL;
        if (!is.null(maf))
          tmp = rbind(tmp, genematrix.mut);
        
        if (!is.null(seg))
          tmp = rbind(tmp, genematrix.cn)
        else if (!is.null(mat))
          tmp = rbind(tmp, mat)
        
        if (!is.null(dranger))
          tmp = rbind(tmp, genematrix.ra);
        
        genes = unique(rownames(tmp)[rowSums(tmp != 0, na.rm = TRUE)!=0])
        gene.ord = intersect(gene.ord, genes)
        patients = colnames(tmp)[colSums(tmp != 0, na.rm = TRUE)!=0]
        
        if (!is.null(maf))
          genematrix.mut = genematrix.mut[intersect(genes, rownames(genematrix.mut)), patients, drop = FALSE]
        if (!is.null(seg))
            genematrix.cn = geneb
        else if (!is.null(mat))
          mat = mat[intersect(genes, rownames(mat)), patients, drop = FALSE]
        
        if (!is.null(dranger))
          genematrix.ra =  genematrix.ra[, patients, drop = FALSE]        
      }

    num.pat = length(patients)
    
    # order patients in "staircase" pattern based on mutation, cn, ra results 
    pat.ord = patients;
    if (!keep.pat.ord)
      {
        tmp.mat = NULL;
        if (!is.null(dranger))
          tmp.mat = rbind(tmp.mat, genematrix.ra)
        
        if (!is.null(maf))
          tmp.mat = rbind(tmp.mat, genematrix.mut[gene.ord,, drop = FALSE]>=MIN.NS) ## order using non-silent coding variants (hacky coding)
        
        if (!is.null(seg))
          tmp.mat = rbind(tmp.mat, genematrix.cn[gene.ord,, drop = FALSE])
        else if (!is.null(mat))
          tmp.mat = rbind(tmp.mat, mat[gene.ord,, drop = FALSE])

        pat.ord = rev(colnames(tmp.mat)[border(t(tmp.mat[(intersect(genes.og, rownames(tmp.mat))), ]!=0))])
      }

    ## listify gene.ord according to gene.lists 
    gene.ord = lapply(gene.lists, function(x) intersect(gene.ord,x))


    #########
    # Prep patient tracks (if exist)
    #########
    if (!is.null(pat.tracks)) # matricize and listify
      {
        if (is.matrix(pat.tracks))
          pat.tracks = as.data.frame(pat.tracks, stringsAsFactors = FALSE)
        
        if (!is.list(pat.tracks) | inherits(pat.tracks, 'data.frame'))
          {
            rnames = rownames(pat.tracks);
            pat.tracks = as.list(pat.tracks)
            tmp = lapply(1:length(pat.tracks), function(x) names(pat.tracks[[x]]) <<- rnames) 
          }
        
        for (i in 1:length(pat.tracks))
          {
            if (is.vector(pat.tracks[[i]]))
              {
                tmp.tracks = as.data.frame(matrix(pat.tracks[[i]], ncol = 1))
                rownames(tmp.tracks) = names(pat.tracks[[i]])
                colnames(tmp.tracks) = names(pat.tracks)[[i]]
                pat.tracks[[i]] = as.matrix(tmp.tracks[pat.ord, , drop = FALSE])
              }
            if (length(intersect(patients, rownames(pat.tracks[[i]])))>0)
              {
                if (is.matrix(pat.tracks[[i]]))
                  pat.tracks[[i]] = t(as.matrix(as.data.frame(pat.tracks[[i]])[pat.ord,, drop = FALSE]))
                else if (is.data.frame(pat.tracks[[i]]))
                  pat.tracks[[i]] = t(as.matrix(pat.tracks[[i]][pat.ord, ,drop = FALSE])) # roundabout way of transposing and creating "NA" entries
                else
                  {
                    warning('Unrecognized patient track format')
                    pat.tracks[[i]] = NULL
                  }
              }
            else # we don't know what it is and flag for removal
              {
                warning('pat.track with empty rownames or no intersection with patient set: will be ignored.')
                pat.tracks[[i]] = NULL;
              }
          }

        pat.tracks = pat.tracks[!sapply(pat.tracks, is.null)]

        if (length(pat.tracks)==0)
            pat.tracks = NULL
      }

    ## if col.pat.tracks provided as single vector, then replicate over all the pat.tracks (will change in the future)
    if (is.character(col.pat.tracks) | c(!is.list(col.pat.tracks) & is.vector(col.pat.tracks)))
      col.pat.tracks = rep(list(col.pat.tracks), length(pat.tracks))
    
    #########
    # BUILD LAYOUT for variant and patient characteristic matrices
    # 
    #########;
    LINES.PER.INCH = 7.575757; ## R constant
    num.rows = sum(sapply(gene.lists, length))
    num.tracks = length(gene.lists) + length(pat.tracks)

    if (is.null(top.marg.lim)) ## set the x limit for the bar plot / marginal axis (e.g. number of patients, or largest number of patients with a given gene or patient attribute 
      {
        top.marg.lim = 1;
        if (!is.null(maf))
          top.marg.lim = max(top.marg.lim, rowSums(genematrix.mut>0))
        else if (!is.null(seg))
          top.marg.lim = max(top.marg.lim, rowSums(genematrix.cn>0))
        else if (!is.null(mat))
          top.marg.lim = max(top.marg.lim, rowSums(!is.na(mat)))
        else if (!is.null(dranger))
          top.marg.lim = max(top.marg.lim, rowSums(genematrix.ra>0))        
      }

    if (is.null(bottom.marg.lim))
      if (!is.null(pat.tracks))
          bottom.marg.lim = max(bottom.marg.lim, unlist(sapply(pat.tracks, function(x) rowSums(!is.na(x)))))

    if (single.scale)
      bottom.marg.lim = top.marg.lim = max(bottom.marg.lim, top.marg.lim)
    
    if (length(pat.tracks)>0) num.rows = num.rows + sum(sapply(pat.tracks, nrow))
    if (!is.null(dranger))
      {
        num.rows = num.rows + nrow(genematrix.ra)
        num.tracks = num.tracks + 1;
      }

    if (pat.labels)
      PAT.LABEL.SIZE = 2
    else
      PAT.LABEL.SIZE = 1;
    
    layout.mat = NULL;
    layout.heights = NULL;
    pc.legend = min(max(0, pc.legend), 100)
    pc.sidebar = min(max(0, pc.sidebar), 100-pc.legend)

    layout.widths = 8*c(1-(pc.sidebar/100)-(pc.legend/100), pc.sidebar/100,  pc.legend/100);
#      layout.widths = 8*c(1-pc.legend/100, pc.legend/100); 
    
    # add a mutation / cn panel for each gene.list (height scaled to number of rows)
    num.mut.tracks = 0
    for (gene.list in gene.lists)    
      if (!is.null(maf) | !is.null(seg) | !is.null(mat))
        {
          layout.mat = rbind(layout.mat, matrix(c(1, 2, NA) + suppressWarnings(sum(!is.na(layout.mat))), nrow = 1, ncol = 3))
          layout.heights = c(layout.heights, pad.track + length(gene.list))
          num.mut.tracks = num.mut.tracks + length(gene.list)
        }
       
    if (length(gene.lists)>0)
      {
        layout.mat[1, 3] = suppressWarnings(sum(!is.na(layout.mat)))+1;
        sum.gene.heights = sum(layout.heights);  # sum of the "gene" Y-real estate on the plot .. helps us position the variant type legends
      }
    
    # if rearrangements exist add rearrangement panel
    num.ra.tracks = 0;
    if (!is.null(dranger))
      {
        layout.mat = rbind(layout.mat, matrix(c(1, 2, 3)+sum(!is.na(layout.mat)), nrow = 1, ncol = 3))
        layout.heights = c(layout.heights, pad.track+nrow(genematrix.ra))
        num.ra.tracks = num.ra.tracks + nrow(genematrix.ra)
      }
    
    # add matrix and legend panel for each patient track (height scaled to number of rows)
    num.pat.tracks = 0;
    for (pat.track in pat.tracks)
      {
        layout.mat = rbind(layout.mat, matrix(c(1, 2, 3)+sum(!is.na(layout.mat)), nrow = 1, ncol = 3))
        layout.heights = c(layout.heights, pad.track+nrow(pat.track))
        num.pat.tracks = num.pat.tracks + nrow(pat.track)
      }

    # add title panel and histogram axis

    # size legend / padding to skinnify mutrix tracks
    # fatness is ratio of height to width
    FATNESS = 0.08;  
    gene.panel.width = (1-(pc.legend/100)-(pc.sidebar/100))*par('din')[1]*(1-(pc.ylab/100))
    tot.tracks = (num.pat.tracks + num.ra.tracks + num.mut.tracks)
    pc.tracks = FATNESS * 100 / ((par('din')[2] / tot.tracks) / gene.panel.width)
    pc.tracks = max(20, min(70, pc.tracks)) # but not too skinny
    
    Y.MAR.GENE = pad.track/2*par('din')[2]*LINES.PER.INCH*(pc.tracks/100)/(num.rows+(num.tracks*pad.track)); # calibrate pad.track and actual mar values we will pass down
    X.LEG = 0.15;
    CEX.X = 0.8; # width of little rectangles (should be less than 1)
    
    top.bar.height = sum(layout.heights)*(100-pc.tracks)/(pc.tracks)
    layout.mat = rbind(matrix(c(1, 2, NA)+sum(!is.na(layout.mat)), nrow = 1, ncol = 3), layout.mat)
    layout.heights = c(top.bar.height*0.5,  layout.heights)
    
    # add label panel on bottom 
    layout.mat = rbind(layout.mat, matrix(c(1, 2, NA)+sum(!is.na(layout.mat)), nrow = 1, ncol = 3))
    layout.heights = c(layout.heights, 0.5*top.bar.height)

    inches.genes = (sum.gene.heights / sum(layout.heights)) * par('din')[2];
    
    # make layout
    layout.mat[is.na(layout.mat)] = 0;
    layout(layout.mat, heights = layout.heights, widths = layout.widths);

    ######
    ## START DRAWING following order of layout matrix
    ######
    pc.ylab = min(max(0, pc.ylab), 100)
    LEFT.PANEL.BOUNDARY = layout.widths[1]/sum(layout.widths)*par('din')[1]*LINES.PER.INCH*(pc.ylab/100);  # translate pc.ylab into R graphics mumbo jumbo
    RIGHT.PANEL.BOUNDARY = 1;
    TOP.BOTTOM.MAR = 2; 
    
    mafseg.legends = list();
    
   # define CN colorway (will be used below)
    if (!is.null(seg))
      {
        cn.var.types = c('Broad Amp', 'Focal Amp', 'Amp + Del', 'Broad Del', 'Focal Del')
        cn.brewer.set = 'Set2';
        cn.var.colors = brewer.pal(8, cn.brewer.set)[c(2, 4, 6, 1, 3)];
        cn.var.colors[3] = 'lightgreen'
        cn.var.colors = lighten(cn.var.colors, 30)
        cn.col.ord = match(cn.map[as.character(1:5)], cn.var.types); # we need this "order" because the order of the legend is not defined in the same order as the numeric indices comprising the cn.matrix (fix?)
        mafseg.legends[['CNA']] = data.frame(label = cn.var.types, border = cn.var.colors, fill = cn.var.colors, stringsAsFactors = F)
      }
    else if (!is.null(mat))      
      if (is.character(mat[1,1])) ## then make standard (CNA) legend
        {
          mafseg.legends[['CNA']] = data.frame(label = names(mat.colormap), border = mat.colormap, fill = mat.colormap, stringsAsFactors = F)
          mat.plot = mat;
        }
      else ## specify spectrum "legend"
        {
          mafseg.legends[['mat']] = data.frame()
          
          ## map spectrum onto mat colormap
          mat.range = range(mat, na.rm = T)
          mat.breaks = seq(mat.range[1]-diff(mat.range)*1e-7, mat.range[2], length.out = length(mat.colormap)+1)

          mat.plot = matrix(cut(mat, mat.breaks, labels = FALSE), nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat)) ## convert to "char" matrix
        }
    
    ##
    ## draw matrix for each genelist
    ##
    
    for (i in 1:length(gene.lists))
      {        
        mafseg.canvas.drawn = FALSE;  ## keep track to see if primary canvas has been drawn for this panel
        par(mar=c(1*Y.MAR.GENE,LEFT.PANEL.BOUNDARY,1*Y.MAR.GENE,RIGHT.PANEL.BOUNDARY), las=1)
        
        if (!is.null(seg))
          {
            this.genematrix.cn = genematrix.cn[gene.ord[[i]], pat.ord, drop = FALSE];
            if (show.plot)
              {
                if (!mafseg.canvas.drawn)
                  {
                    ##
                    ## NEW PANEL
                    ## rectrix
                    recTrix(this.genematrix.cn*0, col = 'gray85', new.plot = T, lab.y = T, col.canvas = col.canvas, col.canvas.border = col.canvas.border , border = T, cex.lab = cex.gene.labels, add.canvas = F, cex.x = CEX.X, cex.y = cex.cn, lwd = lwd.border)
                    mafseg.canvas.drawn = T
                  }           
                recTrix(this.genematrix.cn, col = cn.var.colors[cn.col.ord], new.plot = F, cex.x = CEX.X, cex.y = cex.cn, lwd = lwd.border)
              }            
          }
        else if (!is.null(mat))
          {
            this.mat = mat.plot[gene.ord[[i]], pat.ord, drop = FALSE];
            if (show.plot)
              {
                if (!mafseg.canvas.drawn)
                  {
                    ##
                    ## NEW PANEL
                    ## rectrix
                    recTrix(sign(!is.na(this.mat))*0, col = 'gray85', new.plot = T, lab.y = T, col.canvas = col.canvas, col.canvas.border = col.canvas.border , border = T, cex.lab = cex.gene.labels, add.canvas = F, cex.x = CEX.X, cex.y = cex.cn, lwd = lwd.border)
                    mafseg.canvas.drawn = T
                  }           
                recTrix(this.mat, col = mat.colormap, new.plot = F, cex.x = CEX.X, cex.y = cex.cn, lwd = lwd.border)
              }            
          }
        
        if (!is.null(maf))
          {
                        
            this.genematrix.mut = genematrix.mut[gene.ord[[i]], pat.ord, drop = FALSE]
            if (show.plot)
              {
                mut.col = mut.var.colors[match(mut.map[as.character(1:length(mut.map))], names(mut.var.colors))];
                if (!mafseg.canvas.drawn)
                  {
                    ##
                    ## NEW PANEL
                    ## rectrix
                    recTrix(this.genematrix.mut*0, col = mut.col, new.plot = T, lab.y = T, add.canvas = T, col.canvas = col.canvas, col.canvas.border = col.canvas.border,
                            border = col.canvas, add.grid = T, cex.lab = cex.gene.labels, col.grid = 'white', cex.x = CEX.X, cex.y = cex.cn, lwd = lwd.border)
                    mafseg.canvas.drawn = T
                  }                
                tmp = this.genematrix.mut; tmp[tmp==0] = NA;
                
                recTrix(tmp, col = mut.col, new.plot = F, cex.x = CEX.X*.8, cex.y = cex.mut, add.grid = F, border = col.border, lwd = 0.2)
                ## draw gene list label (if exists)
                if (!is.null(names(gene.lists)))
                  {
                    par(xpd = NA);
                    text(-pos.ylab * ncol(this.genematrix.mut), nrow(this.genematrix.mut)/2 + 0.5,  strwrap2(gsub(gsub.ylab[1], gsub.ylab[2], names(gene.lists)[[i]]), nchar.wrap), adj = adj.ylab, srt = srt.ylab, cex = 1.2 * cex.ylab)
                  }

                if (is.factor(maf.subset$Protein_Change))
                  maf.subset$Protein_Change = as.character(maf.subset$Protein_Change)
                
                # overlay protein change if specified
                if (plot.protein.change)
                  {
                    ## this is a little redundant with code 60 below .. need to clean up
                    maf.subset.coding = maf.subset[nchar(maf.subset$Protein_Change)>0, ]
                    tmp.pchange = aggregate(maf.subset.coding$Protein_Change,
                      by = list(maf.subset.coding$patient_name, maf.subset.coding$Hugo_Symbol),
                      FUN = function(x) paste(unique(x), collapse = ", "))
                    tmp.pchange = tmp.pchange[tmp.pchange[,1] %in% colnames(this.genematrix.mut) & tmp.pchange[,2] %in% rownames(this.genematrix.mut), ]
                    tmp.pchange[,1] = match(tmp.pchange[,1], colnames(this.genematrix.mut))
                    tmp.pchange[,2] = match(tmp.pchange[,2], rownames(this.genematrix.mut))
                    tmp.pchange = tmp.pchange[rowSums(!is.na(tmp.pchange[,1:2]))!=0, ]
                    if (nrow(tmp.pchange)>0)
                        lapply(1:nrow(tmp.pchange), function(x) text(tmp.pchange[x,1], nrow(this.genematrix.mut)-tmp.pchange[x,2]+1, gsub(',', '\n', tmp.pchange[x,3]), col = 'white', srt = srt.protein.change, cex = cex.protein.change))
                }
                
                ##
                ## NEW PANEL
                ## mutation sidebar                
                BAR.SPACE = 1-cex.mut;
                par(mar = c(par('mar')[1], 0, par('mar')[3], 0))
                plot.blank()
                par(usr = c(0, top.marg.lim, 0+BAR.SPACE/2, nrow(tmp)*(1+BAR.SPACE)+BAR.SPACE/2))
                ix = which(!is.na(tmp), arr.ind = T);
                tmp.counts = t(table(rownames(ix[,1, drop = FALSE]), tmp[ix]));
                counts = matrix(0, ncol = nrow(tmp), nrow = length(mut.map)-1, dimnames = list(as.character(1:(length(mut.map)-1)), rev(rownames(tmp))))
                counts[rownames(tmp.counts), colnames(tmp.counts)] = tmp.counts;
                barplot(counts, horiz = TRUE, border = NA, lwd = 0.5,  col = mut.col, axes = FALSE, names.arg = rep("", ncol(counts)), add = T, offset = 0, space = BAR.SPACE)
                lines(c(0, 0), c(par('usr')[3:4]), col = 'black', lwd = 0.5);

                if (!is.null(mark.genes))
                  {
                    if (is.null(names(mark.genes)))
                      {
                        mark.genes = intersect(genes, mark.genes)
                        tmp = rep(8, length(mark.genes));
                        names(tmp) = mark.genes;
                        mark.genes = tmp;
                      }
                    
                    ix = nrow(this.genematrix.mut)-match(names(mark.genes), rownames(this.genematrix.mut))+1;
                    step = (par('usr')[4] - par('usr')[3])/(nrow(this.genematrix.mut))
                    y.coord = seq(par('usr')[3], par('usr')[4], step)[ix]+step/2
                    points(colSums(counts)[ix]+0.07*top.marg.lim, y.coord, pch = mark.genes, cex = cex.gene.labels*0.8)
                  }
                
              }
          }
        else if (show.plot)
          {
            par(mar = c(Y.MAR.GENE, 0, Y.MAR.GENE, 0))
            ##
            ## NEW PANEL
            ## blank sidebar
            plot.blank() 
          }
        
      }
    
    # finalize variant matrices for output (i.e. replace the variant codes with the actual variant type text (or the protein change)
    if (!is.null(seg))
      {
        lapply(names(cn.map), function(x) genematrix.cn <<- array(sub(x, cn.map[x], genematrix.cn), dim = dim(genematrix.cn), dimnames = dimnames(genematrix.cn)))
        genematrix.cn = genematrix.cn[rev(unlist(gene.ord)), pat.ord];
      }
    else if (!is.null(mat))
      mat = mat[rev(unlist(gene.ord)), pat.ord];
        
    if (!is.null(maf))
      {
        tmp = aggregate(maf.subset$Protein_Change, by = list(maf.subset$patient_name, maf.subset$Hugo_Symbol), FUN = function(x) paste(unique(x), collapse = ", "))
        tmp = tmp[tmp[,1] %in% colnames(this.genematrix.mut) & tmp[,2] %in% rownames(this.genematrix.mut), ]

        tmp[,1] = match(tmp[,1], colnames(genematrix.mut))
        tmp[,2] = match(tmp[,2], rownames(genematrix.mut))
        tmp = tmp[rowSums(!is.na(tmp[,1:2]))!=0, ]                
        
        if (dump.protein.change)
          {
            lapply(1:nrow(tmp), function(x) genematrix.mut[tmp[x,2], tmp[x,1]] <<- tmp[x,3])
            genematrix.mut[genematrix.mut=="0"] = "";
          }
        else # dumps only variant types
          lapply(names(mut.map), function(x) genematrix.mut <<- array(sub(x, mut.map[x], genematrix.mut), dim = dim(genematrix.mut), dimnames = dimnames(genematrix.mut)))
        
        mafseg.legends[['Mutation']] = data.frame(label = names(mut.var.colors), border = mut.var.colors, fill = mut.var.colors, stringsAsFactors = F)
        genematrix.mut = genematrix.mut[rev(unlist(gene.ord)), pat.ord, drop = FALSE];
      }
    
    ## draw mut / cn legend
    if (show.plot & length(gene.lists)>0)
      {
        ##
        ## NEW PANEL
        ## mut / seg legend
        
        par(mar=c(0,.001, 0, 0.5))
        par(xpd = NA);
        
        plot.blank(xlim = c(0, 1), ylim = c(-par('pin')[2], 0))
#        else ## hack to allow drawing of a numeric spectrum colormap via rectrix (below)
#          plot.blank(xlim = c(-length(mat.colormap)*0.1, length(mat.colormap)*1.2), ylim = c(-par('pin')[2], 0))
        
        mafseg.legends = rev(mafseg.legends)
         if (length(mafseg.legends)==1)
           {
             leg.loc = list(c(X.LEG*diff(par('usr')[1:2])+ par('usr')[1], -0 * inches.genes))
             leg.yadj = 1
           }
         else if (!is.null(mat))
           {
             leg.loc = list(c(X.LEG*diff(par('usr')[1:2])+ par('usr')[1], -.45*inches.genes), c(X.LEG*diff(par('usr')[1:2])+ par('usr')[1], -0.75*inches.genes))
             leg.yadj = c(0, 1)
           }
         else
           {
             leg.loc = list(c(X.LEG*diff(par('usr')[1:2])+ par('usr')[1], -.45*inches.genes), c(X.LEG*diff(par('usr')[1:2])+ par('usr')[1], -0.55*inches.genes))
             leg.yadj = c(0, 1)
           }
        
        for (i in 1:length(leg.loc))
          if (names(mafseg.legends)[i] == 'mat') ## make numeric spectrum
            {
              Y.PAD = 0.05*inches.genes
              X.WID = 0.2;
#              x = seq(leg.loc[[i]][1]-X.WID, leg.loc[[i]][1]+X.WID, length.out = length(mat.colormap)+1)
              x = seq(leg.loc[[i]][1]+X.WID, leg.loc[[i]][1]+3*X.WID, length.out = length(mat.colormap)+1)              
              rect(x[1:(length(x)-1)], leg.loc[[i]][2] - Y.PAD, x[2:(length(x))], leg.loc[[i]][2] + Y.PAD, col = mat.colormap, border = NA)
              text(x[1] - diff(range(x))*0.1, leg.loc[[i]][2], signif(min(mat, na.rm = T),2) , srt = 0, cex = cex.legend*0.8, adj = c(1, 0.5))
              text(x[length(x)] + diff(range(x))*0.1, leg.loc[[i]][2], signif(max(mat, na.rm = T),2), srt = 0, cex = cex.legend*0.8, adj = c(0, 0.5))
              text(mean(x), leg.loc[[i]][2]+2*Y.PAD, 'Colorscale', srt = 0, cex = cex.legend, adj = c(0.5, 0))
              
            }
          else
            legend(leg.loc[[i]][1], leg.loc[[i]][2], mafseg.legends[[i]]$label, border = mafseg.legends[[i]]$border, fill = mafseg.legends[[i]]$fill, cex=cex.legend, bty = "n", inset=0,
                   title = paste(names(mafseg.legends)[i], 'Types'), ncol = ncol.legend, xjust = 0, yjust = leg.yadj[i])        
      }
 
    if (!is.null(dranger)) ############ UNDER CONSTRUCTION
      {        
        par(mar=c(1*Y.MAR.GENE,LEFT.PANEL.BOUNDARY,1*Y.MAR.GENE,RIGHT.PANEL.BOUNDARY), las=1)
        genematrix.ra = genematrix.ra[rev(border(genematrix.ra[, pat.ord, drop = FALSE]>0)),, drop = FALSE]
        
        if (show.plot)
          {
            ##
            ## NEW PANEL
            ## dranger plot (needs work)
            image(1:length(pat.ord), 1:nrow(genematrix.ra), t(genematrix.ra[, pat.ord, drop = FALSE]), col=c("grey94", brewer.pal(length(ra.uvariants), "Set3")), zlim=c(0,length(ra.uvariants)), xlab="", ylab="", axes=FALSE)
            axis(2, at = 1:nrow(genematrix.ra), labels = rownames(genematrix.ra), tick = FALSE, las = 2,  cex.axis = cex.gene.labels)
            segments(.5 + 1:(length(pat.ord)-1), .5, .5 + 1:(length(pat.ord)-1), .5 + length(gene.ord), col="grey96", lwd=.5)
            segments(.5, .5 + 1:(length(gene.ord)-1), .5 + length(pat.ord), .5 + 1:(length(gene.ord)-1), col="grey96", lwd=.5)
          }
        
        genematrix.ra[genematrix.ra==0] = NA;
        genematrix.ra = array(ra.uvariants[genematrix.ra], dim = dim(genematrix.ra), dimnames = dimnames(genematrix.ra))
        
        if (show.plot)
          {
            ##
            ## NEW PANEL
            ## dranger legend
            par(mar= c(Y.MAR.LEGEND,LEFT.PANEL.BOUNDARY, Y.MAR.LEGEND, 0.5))
            plot(0, type="n", axes=FALSE, xlab="", ylab="");       
            smartlegend("left", "center", ra.uvariants, fill=brewer.pal(length(ra.uvariants), "Set3"), cex=cex.legend, ncol=1, bty = "n", inset=0, title = "Rearrangement Types", pt.cex = 0.75)

            ## sidebar
            par(mar = c(Y.MAR.GENE, 0, Y.MAR.GENE, 0))
            plot.blank()
            text(0, 0, 'Side bar 1')
          }
        
        genematrix.ra = genematrix.ra[, pat.ord, drop = FALSE];
      }

    ##
    ## draw patient track matrices
    ##
    if (!is.null(pat.tracks) & length(pat.tracks)>0)
      {
        for (i in 1:length(pat.tracks))
          {                          
            if (show.plot)
              {
                ##
                ## NEW PANEL
                ## patient track rectangles
                
                par(mar=c(1*Y.MAR.GENE,LEFT.PANEL.BOUNDARY,1*Y.MAR.GENE,RIGHT.PANEL.BOUNDARY), las=1)
                par(xpd = NA)
                uval = setdiff(unique(as.vector(pat.tracks[[i]])), NA)

                # get colorways from col.pat.tracks                
                ucol = col.pat.tracks[[((i-1) %% (length(col.pat.tracks)))+1]]                

                if (is.list(ucol)) ## since it's temping for users to provide named vectors as lists
                  ucol = unlist(ucol)

                if (is.character(uval))
                  {
                    if (!is.null(names(ucol))) ## if ucol has names (ie comes from a custom colorway) , then we assume that pat.tracks[[i]] is a matrix of those names
                      {
                        this.pat.track = pat.tracks[[i]]                        
                        
                        non.overlap = setdiff(setdiff(uval, names(ucol)), NA)
                        if (length(non.overlap)>0) # warn people if they are doing something wrong 
                          warning('There are items in this patient track that do not map to the custom colorway: ',
                                  paste(non.overlap, collapse = ","))

                        uval = names(ucol)
                        this.pat.track[this.pat.track %in% non.overlap] = NA;
                      }
                    else # ucol is unnamed, so we just map uval[k] to ucol[k]
                      {
                        ucol = ucol[(1:length(uval) %% length(ucol))+1]
                        this.pat.track = matrix(match(pat.tracks[[i]], uval), nrow = nrow(pat.tracks[[i]]), ncol = ncol(pat.tracks[[i]]))
                      }
                    
                    dimnames(this.pat.track) = dimnames(pat.tracks[[i]])
                    rownames(this.pat.track) = gsub(gsub.ylab[1], gsub.ylab[2], rownames(this.pat.track)); # remove unwanted characters from labels
                    recTrix(this.pat.track, col = ucol, border = col.border, cex.x = CEX.X, cex.y = cex.cn, lwd = lwd.border, new.plot = T, lab.y = T, add.canvas = F)
                  }
                else # this pat track is numeric or logical 
                  {
                    this.pat.track = pat.tracks[[i]];
                    dimnames(this.pat.track) = dimnames(pat.tracks[[i]])                    
                    ucol = ucol[1];
                    if (is.numeric(uval))
                      recTrix(this.pat.track, col = ucol, border = col.border, cex.x = CEX.X, cex.y = cex.cn, lwd = lwd.border, new.plot = T, lab.y = T, add.canvas = F)
                    else if (is.logical(uval))
                      recTrix(this.pat.track, col = col.true, border = col.border, cex.x = CEX.X, cex.y = cex.cn, lwd = lwd.border, new.plot = T, lab.y = T, add.canvas = F)
                    else
                      plot.blank(xlim = c(0, 1), ylim = c(-1, 1))
                  }

               # draw label (if exists)
                if (!is.null(names(pat.track)) & nrow(pat.track)>1)
                  {
                    par(xpd = NA);
                    text(-pos.ylab * ncol(pat.tracks[[i]]), nrow(pat.tracks[[i]])/2 + 0.5,  strwrap2(gsub(gsub.ylab[1], gsub.ylab[2], names(pat.tracks)[[i]]), nchar.wrap), adj = c(0.5, 0), srt = 90, cex = 1.2 * cex.ylab)
                  }
                
                par(xpd = NA);
                
                if (length(uval)>0)
                  {
                    ## sidebar
                    BAR.SPACE = 0.2;
                    
                    ##
                    ## NEW PANEL
                    ## patient track legend
                    
                    par(mar = c(par('mar')[1], 0, par('mar')[3], 0))
                    plot.blank()
                    par(usr = c(0, bottom.marg.lim, 0+BAR.SPACE/2, nrow(this.pat.track)*(1+BAR.SPACE)+BAR.SPACE/2))

                    if (is.character(uval))  ## create stacked barplot of unique colors 
                      {
                        if (!is.null(names(ucol)))
                          counts = matrix(0, ncol = nrow(this.pat.track), nrow = length(uval), dimnames = list(uval, rev(rownames(this.pat.track))))
                        else
                          counts = matrix(0, ncol = nrow(this.pat.track), nrow = length(uval), dimnames = list(as.character(1:length(uval)), rev(rownames(this.pat.track))))
                                                
                        ix = which(!is.na(this.pat.track), arr.ind = T);
                        tmp.counts = t(table(names(ix[,1]), this.pat.track[ix]))
                        counts[rownames(tmp.counts), colnames(tmp.counts)] = tmp.counts;
                        barplot(counts, horiz = TRUE, border = NA, lwd = 0.5,  col = ucol, axes = FALSE, names.arg = rep("", ncol(counts)), add = T, offset = 0, space = BAR.SPACE)
                      }
                    else if (is.numeric(uval)) ## create stacked "spectrum" barplot summarizing numeric distribution
                      {
                        order.ix = order(this.pat.track[1,]);
                        this.pat.track.colors = col.scale(this.pat.track[1,order.ix], val.range = range(this.pat.track[1, !is.na(this.pat.track[1,])]), col.max = ucol)
                        scaled.ucol = unique(this.pat.track.colors);
                        counts = matrix(0, ncol = nrow(this.pat.track), nrow = length(scaled.ucol), dimnames = list(scaled.ucol, rev(rownames(this.pat.track)))) # counts is scaled.ucol x individuals
                        ix = which(!is.na(this.pat.track[1,order.ix, drop = FALSE]), arr.ind = T);
                        tmp.counts = t(table(names(ix[,1]), this.pat.track.colors[ix[,2]]))                        
                        counts[rownames(tmp.counts), colnames(tmp.counts)] = tmp.counts;
                        barplot(counts, horiz = TRUE, border = NA, lwd = 0.5,  col = unique(this.pat.track.colors), axes = FALSE, names.arg = rep("", ncol(counts)), add = T, offset = 0, space = BAR.SPACE)
                      }
                        
                    lines(c(0, 0), c(par('usr')[3:4]), col = 'black', lwd = 0.5);
                   
                    if (is.character(uval))
                      {
                        par(mar = c(par('mar')[1], 0.001, par('mar')[3], 0.5))
                        plot.blank(xlim = c(0, 1), ylim = c(0, 1))
                        legend(X.LEG, 0.5, uval, border = ucol, fill= ucol, pt.cex=cex.legend*0.8, cex = cex.legend, bty = "n", inset=0,
                               ncol = max(2, ceiling(length(uval)/2)), xjust = 0, yjust = 0.5)
                      }
                    else if (is.numeric(uval) & (any(!is.na(uval))))
                      {
                        par(mar = c(par('mar')[1], 0.001, par('mar')[3], 0.5))
                        NUM.SHADES = 10
                        plot.blank(xlim = c(-X.LEG*(NUM.SHADES+5)*2,NUM.SHADES+5), ylim = c(-1, 3))
                        col.scale = matrix(seq(min(uval), max(uval), diff(range(uval))/NUM.SHADES), nrow = 1)
                        recTrix(col.scale, col = ucol[1])
                        text(0.2, 1, range(uval)[1], srt = 0, cex = cex.legend*0.8, adj = c(1, 0.5))
                        text(length(col.scale)+0.8, 1, range(uval)[2], srt = 0, cex = cex.legend*0.8, adj = c(0, 0.5))
                      }
                    else if (is.logical(uval))
                      {
                        par(mar = c(par('mar')[1], 0.001, par('mar')[3], 0.5))                    
                        plot.blank(xlim = c(0, 1), ylim = c(0, 1))
                        legend(X.LEG, 0.5, string.true, border = col.true, fill = col.true, pt.cex=cex.legend*0.8, cex = cex.legend, bty = "n", inset=0,
                               ncol = max(2, ceiling(length(uval)/2)), xjust = 0, yjust = 0.5)
                      }
                    else
                      {
                        plot.blank(xlim = c(0, 1), ylim = c(0, 1))
                      }
                  }
                else
                  {
                    par(mar = c(0,0,0,0));
                    plot.blank();
                    plot.blank();
                  }
              }        
          }
      }

    ## draw top bar plot (if exists)
    if (show.plot)
      {        
        if (!is.null(bar.top))
          {
            ##
            ## NEW PANEL
            ## bar plot on top
            
            par(mar=c(0,LEFT.PANEL.BOUNDARY,2*TOP.BOTTOM.MAR,RIGHT.PANEL.BOUNDARY), las=1, xpd = TRUE)
            plot.blank()            
            if (!is.null(bar.top$data))
              {
                BAR.SPACE = 0.2;

                if (ncol(bar.top$data)==1) ## convenience function if we provide a single column of data --> we assume data is transposed
                   bar.top$data = t(bar.top$data)
                
                bar.top$data = as.matrix(bar.top$data)
                bar.top$height = matrix(0, nrow = nrow(bar.top$data), ncol = length(pat.ord), dimnames = list(rownames(data), pat.ord));
                col.ix = intersect(pat.ord, colnames(bar.top$data))
                bar.top$height[, col.ix] = bar.top$data[, col.ix]

                if (is.null(bar.top$log))
                  plot.log = FALSE
                else
                  if (bar.top$log %in% c('y', 'xy', 'yx'))
                    plot.log = TRUE
                  else
                    plot.log = FALSE
                
                bar.top$log = NULL;

                if (plot.log)
                  bar.top$height = log10(bar.top$height)
                
                if (is.null(bar.top$ylim))
                  if (plot.log)
                    bar.top$ylim = c(0, max(bar.top$height)*1.2)
                  else
                    bar.top$ylim = c(0, max(bar.top$height)*1.2)
                else
                  if (plot.log)
                    {
                      if (any(bar.top$ylim<0))
                        warning('Negative values cannot be plotted on a log scale')
                      bar.top$ylim = log10(bar.top$ylim)
                    }
                
                par(usr = c(0 + BAR.SPACE/2, length(pat.ord)*(1+BAR.SPACE) + BAR.SPACE/2,  bar.top$ylim))                
                bar.top$ylim = bar.top$data = NULL;
                if (is.null(bar.top$col))
                  if (nrow(bar.top$height)==1)
                    bar.top$col = 'gray90'
                  else
                    bar.top$col = brewer.pal(8, 'Dark2')[(1:nrow(bar.top$height) %% 8)+1]
                bar.top$add = TRUE;
                bar.top$axes = FALSE;

                if (is.null(bar.top$cex.main))
                  bar.top$cex.main = 0.8;
                bar.top$names.arg = rep('', ncol(bar.top$height))
                bar.top$space = BAR.SPACE;

                at = bar.top$at;
                bar.top$at = NULL;
                
                do.call('barplot', bar.top)
                
                if (!plot.log)
                  {
                    if(is.null(at))                      
                      at = pretty(bar.top$height, 3)
                                        
                    axis(2, at = at, cex.axis = cex.tick*0.8, line = 0.5)
                  }
                else 
                  {
                    if (is.null(at))
                      {
                        at = pretty(10^(bar.top$height), 3);
                        at = seq(round(min(bar.top$height)), max(bar.top$height), 1);
                      }
                    else
                      at = log10(at);                   

                    axis(2, at = at, labels = 10^(at),  cex.axis = cex.tick*0.8, line = 0.5)
                  }                

              }
            else
                warning('No data field in bar.top, ignoring')
          }
        else
          plot.blank()
      }      
      
    ##
    ## NEW PANEL
    ## top sidebar axis
    if (show.plot)
      {
        par(mar = c(0, 0, TOP.BOTTOM.MAR, 0))
        plot.blank();
        
        if (!is.null(maf))
          {
            if (freq.sidebar)
              {
                leg.text = '% prevalence'
                top.marg.lim = top.marg.lim/num.pat*100
              }
            else
              leg.text = '# Patients'
            
            par(usr = c(0, top.marg.lim, 0, 1), xpd = NA)
            axis(3, pos = 0, cex.axis = cex.tick*0.8)
            
            text(mean(par('usr')[1:2]), min(1, par('pin')[1]*0.6/par('pin')[2]), leg.text, adj = c(0.5, 1), cex = cex.tick*0.8)
          }
        
        
        ##
        ## NEW PANEL
        ## bottom label(s)     
        
        if (pat.labels)
          {
            par(mar=c(TOP.BOTTOM.MAR,LEFT.PANEL.BOUNDARY,0,RIGHT.PANEL.BOUNDARY), las=1)
            plot.blank(xlim = c(0.5, length(pat.ord)+0.5), ylim = c(-1, 1))
            adj = c(0, 0.5)
            text((1:length(pat.ord)), rep(1-adj[1], length(pat.ord)), pat.ord, srt = srt.pat.labels, cex = cex.pat.labels, adj = adj)
          }
        else
          {
            par(mar=c(TOP.BOTTOM.MAR,LEFT.PANEL.BOUNDARY,0,RIGHT.PANEL.BOUNDARY), las=1)
            plot.blank(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5))          
          }
        
        ##
        ## NEW PANEL
        ## bottom sidebar axis        
        par(mar = c(TOP.BOTTOM.MAR, 0, 0, 0))    
        if (!is.null(bottom.marg.lim) & !is.null(maf))
          {
            if (freq.sidebar)
              {
                leg.text = '% prevalence'
                bottom.marg.lim = bottom.marg.lim/num.pat*100
              }
            else
              leg.text = '# Patients'
            
            plot.blank(xlim = c(0, bottom.marg.lim), ylim = c(0, 1))
            par(usr = c(0, bottom.marg.lim, 0, 1), xpd = NA)
            axis(1, pos = 1, cex.axis = 0.8*cex.tick)
            text(mean(par('usr')[1:2]), max(0, 1-par('pin')[1]*0.6/par('pin')[2]), leg.text, adj = c(0.5, 0), cex = 0.8*cex.tick)
          }
        else
          plot.blank()
      }
    
    ### put together output matrix 
    out = list();
    if (!is.null(maf))
      out[['mut']] = genematrix.mut;
    if (!is.null(seg))
      out[['cn']] = genematrix.cn;
    if (!is.null(mat))
      out[['mat']] = mat;
    if (!is.null(dranger))
      out[['dranger']] = genematrix.ra;
    if (length(out)==1)
      out = out[[1]];
    
    return(out)    
  }


#################################
# quickTrix
#
# Quick nicely formatted matrix of tabular data in mat (untransposed), with borders drawn around each line
# "..." get passed to "image"
#
#################################
quickTrix = function(mat, label.ticks = 'xy', # if 'x' only draws x-ticks, 'y' only draws y ticks 
  xlab = "", ylab = "", cex.lab = 1, cex.lab.x = cex.lab, cex.lab.y = cex.lab, 
  lwd.border = 1, x.side = 'bottom', y.side = "left", pad = 8,
  mar = NULL,
  y.adj = 0, # adjustment for row labels
  x.adj = 0, # adjustment for col labels
  asp = 1,  # 1 = square, < 1 fat, > 1 skinny
  ...)
  {
    if (is.null(mar))
      par(mar= pad*rep(1, 4))
    else
      par(mar = mar)

    asp = ncol(mat)/nrow(mat)*asp;
    
    image(1:ncol(mat), 1:nrow(mat), t(mat), xlab = xlab, ylab = ylab, axes = F, asp = asp, ...)
    segments(.5 + 1:(ncol(mat)-1), .5, .5 + 1:(ncol(mat)-1), .5 + nrow(mat), col="grey96", lwd=lwd.border)
    segments(.5, .5 + 1:(nrow(mat)-1), .5 + ncol(mat), .5 + 1:(nrow(mat)-1), col="grey96", lwd=lwd.border)            

    lab.ex = 0.02;
    if (label.ticks %in% c('x', 'xy'))        
      if (x.side == "top")
        axis(3, at = 1:ncol(mat), labels = colnames(mat), tick = FALSE, las = 2,  cex.axis = cex.lab.x, pos = nrow(mat)+x.adj)
      else
        axis(1, at = 1:ncol(mat), labels = colnames(mat), tick = FALSE, las = 2,  cex.axis = cex.lab.x, pos = x.adj)

    if (label.ticks %in% c('y', 'xy'))        
      if (y.side == "left")
        axis(2, at = 1:nrow(mat), labels = rownames(mat), tick = FALSE, las = 1,  cex.axis = cex.lab.y, pos = y.adj)
      else
        axis(4, at = 1:nrow(mat), labels = rownames(mat), tick = FALSE, las = 1,  cex.axis = cex.lab.y, pos = ncol(mat) + y.adj)
                
  }


##########################
# cn_matrix
#
# Makes matrix of genes vs individuals where each entry aggregates min or max copy number value in a given gene 
#
# Takes in seg data frame which is data frame with columns $ID, $chrom, $loc.start, $loc.end or $ID, $chr, $start, $end, or $ID, $chr, $pos1, $pos2 + $seg.mean
#
# Can optionally specify which patients to include in matrix. can also specify alternate segment value fields to aggregate (ie other than seg.mean)
##########################
cn_matrix = function(genes, seg, patients = NULL, rg = NULL, hg19 = F,
  field = 'seg.mean', # can be any other field in the seg file
  agg = "max", # can also be "min", "mean"
  use.na = F # will put NA for blank regions, if = F then will use 0
  )
{
  seg = standardize_segs(seg)
 
  if (is.null(rg))
    rg = read_refGene(hg19 = hg19);
  
  ix = match(genes, rg$gene_sym);
  genes.df = data.frame(gene = genes, stringsAsFactors = F);  
  genes.df[, c('chr', 'pos1', 'pos2')] = rg[ix, c('chr', 's1', 'e1')]
  rownames(genes.df) = genes.df$gene;  
  these.segs = seg[seg.on.seg(seg, genes.df),];
  
  uid = unique(seg$ID);
  if (use.na)
    mut.mat = array(NA, dim = c(nrow(genes.df), length(uid)),  dimnames = list(genes.df$gene, uid))
  else
    mut.mat = array(0, dim = c(nrow(genes.df), length(uid)),  dimnames = list(genes.df$gene, uid))

  ## collect min and max raw seg value for a segment overlapping a given known gene to call "amp" and "del" events (poor mans copy annotation)
  for (gene in genes.df$gene)
    {
      tmp.seg = these.segs[seg.on.seg(these.segs, genes.df[gene,]),]
      if (nrow(tmp.seg)>0)
        {
          if (agg == "max")
            {
              tmp = aggregate(tmp.seg[, field], list(tmp.seg$ID), function(x) if(any(!is.na(x))) max(x, na.rm = TRUE) else NA);
              mut.mat[gene, tmp$Group.1] = tmp$x
            }
          else if (agg == "min")
            {
              tmp = aggregate(tmp.seg[,field], list(tmp.seg$ID), function(x) if(any(!is.na(x))) min(x, na.rm = TRUE) else NA);
              mut.mat[gene, tmp$Group.1] = tmp$x;
            }
        }
    }

  if (!is.null(patients))
    {
      others = setdiff(patients, colnames(mut.mat))
      mut.mat = cbind(mut.mat, matrix(NA, nrow = nrow(mut.mat), ncol = length(others), dimnames = list(rownames(mut.mat), others)))
      mut.mat = mut.mat[, patients, drop = FALSE];
    }

  return(mut.mat)
}

################################
# recTrix
#
# Primitive to nicely visualize matrix on existing or new plot using rect (instead of image)
# with (scaled) intensities corresponding to provided color (if col is length 1)
# or colormap to unique non-NA integer values of A (if col is vector)
#
# add.canvas will add (default grey) canvas and add.grid will add nico-style grid (default color another shade of gray)
################################
recTrix = function(A, # A is a matrix of integers or reals
  col = 'black', # col can be (1) a single color, in which case the values of A are scaled as shades of that color (where white = lower range) or
                 # (2) a vector of colors, in which A is assumed to contain positive integers that index col, A entries that map to 0 are given a transparent color in this situation
                 # (3) a NAMED vector of colors, in which A is assumed to be composed of NA and names(col)
  fill = TRUE, border = TRUE, # these specify whether the rectangles should be filled or have borders 
  val.range = NULL, # in the case of shading, the val range can be provided to map colors to a prespecified range
  cex = 1, cex.x = cex, cex.y = cex, lwd = 0.5, new.plot = F, cex.lab = 1, lab.x = NULL, lab.y = NULL, pad = 0,
  shift = 0, shift.x = shift, shift.y = shift, # shift will move entire graph in x and / or y direction
  mat.ord = T, # "matrix orientation" if true then lower "i" indices are at the top (ie the way we imagine matrices to look)
  add.canvas = new.plot, # will add canvas in background of matrix plot with color col.canvas
  col.canvas = 'gray91', #  
  col.canvas.border = col.canvas, #
  add.grid = F, # add grid on top of plotted matrix
  col.grid = 'gray99', # color of grid lines
  lwd.grid = 0.7
  )
{
  
  if (new.plot)
    {
      x.pad = pad*ncol(A)
      y.pad = pad*nrow(A)
      xlim = c(1-0.5-x.pad, ncol(A)+0.5+x.pad)
      ylim = c(1-0.5-y.pad, nrow(A)+0.5+y.pad);
      plot.blank(ylim = ylim, xlim = xlim);
      par(usr = c(par('usr')[1:2], ylim))
      par(xpd = NA);
      if (!is.null(lab.y) & !is.null(rownames(A)))
        {
          if (!(lab.y %in% c(2, 4)))
            lab.y = 2;
          axis(lab.y, at = rev(1:nrow(A)), labels = rownames(A), tick = FALSE, las = 2,  cex.axis = cex.lab, line = 0)
        }
          
      if (!is.null(lab.x) & !is.null(colnames(A)))
        {
          if (!(lab.x %in% c(1, 3)))
            lab.x = 1;
          axis(lab.x,  at = 1:ncol(A), labels = colnames(A), tick = FALSE, las = 2, cex.axis = cex.lab, line = 0)
        }      
    }

  if (is.list(col))  ## since it's tempting to provide named vectors as lists 
    col = unlist(col)

  if (new.plot)
    add.canvas = T

  if (add.canvas)
      rect(1-0.5, 1-0.5, ncol(A)+0.5, nrow(A)+0.5, border = col.canvas.border, col = col.canvas)
        
  if (is.logical(A))
    A[!A] = NA;

  if (mat.ord)
    A = A[rev(1:nrow(A)), , drop = FALSE]

  if (sum(!is.na(A))==0)
    return()
  
  if (length(col)==1) # col is length 1 --> assume color scaling of real values
    {
      if (is.null(val.range))
        val.range = range(A[!is.na(A)])
      
      if (diff(val.range)==0)
        col[1:length(which(!is.na(A)))] = col
      else
        col = col.scale(A[!is.na(A)], val.range = val.range, col.max = col, na.col = 'white')
    }
  else # col is vector --> matrix values are interpreted to index col vector
    {
      vals = A[!is.na(A)]
      og.col = col;      
      col = rep(NA, length(vals))
      col[vals!=0] = og.col[vals[vals!=0]]
    }
  
  if (fill)
    fill = col
  else
    fill = NULL;

  if (!is.null(border))
    {
      if (is.logical(border))
      {
        if (border)
          border = col
        else
          border = col.canvas
      }
      else
        {
                                        # border = border[1]; # this interprets border as a scalar color
          vals = A[!is.na(A)]
          og.border = border;      
          border = rep(NA, length(vals))
          border[vals!=0] = og.border[vals[vals!=0]]
        }
    }
  else
    border = col
  
  ix = which(!is.na(A), arr.ind = TRUE);
  
  if (nrow(ix)>0)
    {
      T = 0.5;
      M = t(sapply(1:nrow(ix), function(x) ix[x,] + c(-T*cex.y, -T*cex.x, T*cex.y, T*cex.x)));
      rect(M[,2] + shift.x, M[,1] + shift.y, M[,4] + shift.x, M[,3] + shift.y, col = fill, border = border, lwd = lwd)
    }

  if (add.grid)
    {
      segments(.5 + shift.x + 1:(ncol(A)-1), .5 + shift.y, .5 + 1:(ncol(A)-1) + shift.x, .5 + nrow(A) + shift.y, col= col.grid, lwd=lwd.grid)
      segments(.5 + shift.x, .5 + 1:(nrow(A)-1) + shift.y , .5 + ncol(A) + shift.x, .5 + 1:(nrow(A)-1) + shift.y, col= col.grid, lwd= lwd.grid)           
    }
  
}


#####################################
# spy
#
# utility function visualizes sparsity structure of matrix
# using rectrix call
######################################
spy = function(A, signed = FALSE, col = NULL, col.grid = NA, lwd.grid = 1.2, border = NULL,  ...)
  {
    if (is(A, 'sparseMatrix'))
      A = as(A, 'matrix')
    
    if (is.null(rownames(A)))
      rownames(A) = as.character(1:nrow(A))

    if (is.null(colnames(A)))
      colnames(A) = as.character(1:ncol(A))

    if (signed)
      {
        if (is.null(col))
          col = c('blue', 'red')
                
        if (is.null(names(col)))
          names(col) = c('-1', '1');

        A[,] = as.character(sign(A))
        recTrix(A, new.plot = T, col = col, col.grid = col.grid, add.grid = TRUE, lab.x = 1, lab.y = 2, border = border, lwd.grid = lwd.grid, ...)
      }
    else
      {
        if (is.null(col))
          col = 'black'
        recTrix(A != 0, new.plot = T, col = col, col.grid = col.grid, add.grid = TRUE, lab.x = 1, lab.y = 2, border = border, lwd.grid = lwd.grid,  ...)
      }
    
  }


###########################
# axis.subplot
#
# function to draw simple horizontal and vertical axes at arbitrary positions on plot with
# decoupled data and plot coordinates (useful when there are many subplots within a single layout item)
#
###########################
axis.subplot = function(side = 1, # can be 1 (horizontal) or 2 (vertical)
  data.range, # (data coordinates) min and max data coordinate value for ends of axis line
              # corresponding to ypos (or xpos) plot coordinates in case of vertical (or horizontal) plot, respectively
  n.ticks = 5, # optimal number of ticks using "pretty" function
  at = NULL, # (data coordinates) tick locations same as in original axis function, over-rides n.ticks 
  ypos = NULL, # (plot coordinates) single value for horizontal, two values (range) for vertical correponding to ypos of axis line
  xpos = NULL, # (plot coordinates) single value for vertical, two values (range) for horizontal corresponding to xpos of axis line
  lwd = 1,  
  lwd.ticks = lwd,
  col = 'black',
  col.ticks = col,
  tick.length = 0.1, # in plot coordinates
  tick.orient = 1, # if >0 then facing "positive" direction, otherwise "negative"
  srt.tick = 0, # rotation of tick.labels relative to "default" (which is 0 for horizontal and 90 for vertical axis)
  cex.tick = 1, # size of tick labels
  adj.tick = c(0.5, 0.5), # justification of tick.labels (overrides defaults)
  ...) # arguments to line)
{  
  if (side == 1 & !(length(xpos)==2 & length(ypos)==1))
    stop('ypos must be of length 1 and xpos of length 2, for side = 1')
  else if (side == 2 & !(length(xpos)==1 & length(ypos)==2))
    stop('ypos must be of length 2 and xpos of length 1, for side = 2')
  else if (!(side %in% c(1, 2)))
    stop('side must be = 1 or = 2')
  
  axis.lim = data.frame(x = xpos, y = ypos);

  if (is.null(at))
    at = pretty(data.range, n = n.ticks)

  # draw axis line
  lines(axis.lim$x, axis.lim$y, lwd = lwd, col = col, ...)

  # draw ticks
  if (tick.length <= 0)
    tick.length = -tick.length;
  
  if (side == 1)
    {
      at.plot = coord.subplot(at, data.range = data.range, plot.range = xpos)
      segments(at.plot, rep(ypos, length(at)), y1 = rep(ypos, length(at))+tick.length, lwd = lwd, col = col, ...)
      text(rep(at.plot, 2), rep(ypos, length(at))+1.5*tick.length, at, adj = adj.tick, cex = cex.tick, srt = 0+srt.tick)
    }
  else if (side == 2)
    {      
      at.plot = coord.subplot(at, data.range = data.range, plot.range = ypos)
      segments(rep(xpos, length(at)), at.plot, x1 = rep(xpos, length(at))+tick.length, lwd = lwd, col = col, ...)
      text(rep(xpos, length(at))+1.5*tick.length, rep(at.plot, 2), at, adj = adj.tick, cex = cex.tick, srt = -90+srt.tick)      
    }      
}

##############################
# coord.subplot
#
# Scales data coordinates in given "data.range" to plot coordinates in "plot.range".  Can be used for x or y coordinates.
# Data coordinates that map outside of the data.range are given NA values. 
#
# Useful for doing many subplots on a single plot axis. 
# e.g. if we want to plot genomic positions 54,100,000 to 56,200,000 on plot coordinates 0.1 to 0.4
##############################
coord.subplot = function(data, # any length data vector  
              data.range, # length 2 vector 
              plot.range # length 2 vector
              )
{
  data[data<data.range[1] | data>data.range[2]] = NA;
  return(((data-data.range[1])/diff(data.range))*diff(plot.range)+plot.range[1])
}

############################
# strwrap2
#
# Wrapper around strwrap that returns a vector of characters the same length as input with new lines
# inserted in word breaks to obey a provided string width.
#
#############################
strwrap2 = function(x, width, sep = NULL, newline = '\n', ...)
  {
    return(sapply(strwrap(x, width = width, simplify = F),
                  function(x) paste(x, collapse = newline)))
  }


#################
#' @name lighten
#'
#' lightens / darkens colors by brighness factor f (in -255 .. 255) that will make lighter if > 0 and darker < 0
#'
#' @param col vector of colors
#' @param f numeric brightness factor
#' @return colors lightened or darkened by factor
################
lighten = function(col, f)
  {
    M = col2rgb(col)
    return(apply(matrix(pmax(0, pmin(255, M + f*matrix(rep(1, length(M)), nrow = nrow(M)))), ncol = length(col))/255, 2, function(x) rgb(x[1], x[2], x[3])))
  }



##########################
# plot.blank
#
# Shortcut for making blank plot with no axes
##########################
plot.blank = function(xlim = c(0, 1), ylim = c(0,1), xlab = "", ylab = "", axes = F, ...)
  {
    plot(0, type = "n", axes = axes, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)
#    par(usr = c(xlim, ylim))
  }

border = function(B, na.rm = TRUE)
  {
    B = array(as.logical(B), dim = dim(B))
    tmp = vector(mode = "numeric", length = nrow(B));
    if (na.rm)
      B[is.na(B)] = FALSE;
    for (i in 1:ncol(B))
        tmp = tmp + 2^(ncol(B)-i)*as.numeric(B[,i]==1);
    return(order(tmp))
  }
