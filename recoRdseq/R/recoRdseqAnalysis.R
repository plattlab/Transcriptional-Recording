## Load Packages ##
x<-c("ggplot2", "VennDiagram",  "pheatmap", "ggfortify", "cluster", "DESeq2", "edgeR", "baySeq", "readxl", "reshape2",  "gplots", "RColorBrewer")
invisible(suppressMessages(lapply(x, require, character.only = TRUE)))


recoRdSeqAnalysis  <- function(
  countsMatrix,  # path to featurecounts output
  designMatrix, #  path to designMatrix
  outPath, # path for output
  totalCountsFile = NULL,
  designFormula =  NA,
  nGenes = 50, # no of top genes for gene box plots (counts) and PCA (variance)
  K = 2, # initialization value for K-means/ FANNY clustering
  minCountsPerSample = 500, # minimum number of total counts in a sample or it is excluded from analysis
  geneBoxAndWhiskerPlots = TRUE,
  transformation = "rlog", # Data transformation prior to plotting gene count boxplots and PCA - can be "log2", "rlog", "vst", "tmm" or NULL
  clustering = FALSE,
  PCAplots = TRUE,
  vennDiagrams = TRUE,
  processInput = TRUE
)
{
  data <- as.data.frame(read.table(countsMatrix, header = TRUE))
  for(i in 7:dim(data)[2]){
    colnames(data)[i]<-strsplit(colnames(data)[i], "/")[[1]][length(strsplit(colnames(data)[i], "/")[[1]])]
    colnames(data)[i]<-gsub(".bam", "", colnames(data)[i])
  }
  rownames(data)<-data[,1]
  data<-data[,-c(1:6)]
  data<-data[which(rowSums(data)>as.numeric(quantile(rowSums(data))[2])),]
  no=nGenes

  # Read in design matrix
  design <- as.data.frame(read_excel(designMatrix))
  rownames(design) <- design[,1]
  design <- design[,-1, drop=FALSE]

  if(!dir.exists(outPath)){dir.create(outPath)}

  # Read in total counts by sample
  if(!is.null(totalCountsFile)){
    totalCounts <- as.data.frame(read.table(totalCountsFile, header = TRUE))
    rownames(totalCounts)<-totalCounts[,1]
    totalCounts <- totalCounts[, -1]
    if(processInput) {
      DEList<-.InputProcess(data=data, design=design, totalCounts=totalCounts, minCountsPerSample=minCountsPerSample)
    } else {
      DEList<-list(data, design, totalCounts)
    }
  } else {
    if(processInput) {
      DEList<-.InputProcess(data=data, design=design, totalCounts=NULL, minCountsPerSample=minCountsPerSample)
    } else {
      DEList<-list(data, design)
    }
  }
  data<-as.data.frame(DEList[[1]])
  design<-as.data.frame(DEList[[2]])
  # DEList<-.removeOutliers(data, design)
  if(!is.null(totalCountsFile)){
    totalCounts<-as.data.frame(DEList[[3]])
  }

  .DEStats(DEList, outPath, designFormula, K=K, geneBoxAndWhiskerPlots=geneBoxAndWhiskerPlots, totalCountsFile, PCAplots=PCAplots, clustering=clustering, vennDiagrams=vennDiagrams, no=no, transformation = transformation)
    
}

## UTILITIES

.InputProcess <- function(data, design, totalCounts, minCountsPerSample=100) {
  # Filtering and sorting counts data by samples in design matrix. This assumes the samples are named consistently in all the inputs
  data<-data[,which(colnames(data)%in%rownames(design))]
  data<-data[,match(rownames(design),colnames(data))]
  if(!is.null(totalCounts)){
    totalCounts<-totalCounts[which(rownames(totalCounts)%in%rownames(design)),]
    totalCounts<-totalCounts[match(rownames(design),rownames(totalCounts)),]
    idx <- which(colSums(data)<minCountsPerSample)
    if (length(idx)>0) {
      design<-design[-idx, ,drop=FALSE]
      totalCounts<-totalCounts[-idx,]
      data<-data[,-idx]
    }

    DEList<-list(data, design, totalCounts)
    DEList
  } else {
    idx <- which(colSums(data)<minCountsPerSample)
    if (length(idx)>0) {
      design<-design[-idx, ,drop=FALSE]
      data<-data[,-idx]
    }

    DEList<-list(data, design)
    DEList
  }
}

# .removeOutliers<-function(data, design, Z_max=2.5)
# {
#   design_o<-design[order(design[,1]),]
#   replicates<-c()
#   rep=1
#   for(rw in 1:dim(design_o)[1]){
#     if(rw>1) {
#       if(isTRUE(all.equal(design_o[rw,],design_o[rw-1,],check.attributes = FALSE)))
#       {
#         replicates<-c(replicates,rep)
#       } else {
#         rep=rep+1
#         replicates<-c(replicates,rep)
#       }
#     } else {
#       replicates<-c(replicates,rep)
#     }
#   }
#   data_o<-data[, ,match(rownames(design_o),colnames(data))]
#   rep_rm<-data.frame(replicates, colSums(data_o))
#   rownames(rep_rm)<-rownames(design_o)
#   replicates<-unique(replicates)
#   Z<-c()
#   for(i in 1:length(replicates)){
#     xm<-median(rep_rm[which(rep_rm[,1]==replicates[i]),2])
#     mad<-median(abs(rep_rm[which(rep_rm[,1]==replicates[i]),2]-xm))
#     Z<-c(Z,0.675*(rep_rm[which(rep_rm[,1]==replicates[i]),2]-xm)/mad)
#   }
#   design_o<-design_o[-which(abs(Z)>Z_max),]
#   design<-design[which(rownames(design)%in%rownames(design_o)),]
#   data<-data[,which(colnames(data)%in%rownames(design))]
#   list(data,design)
# }
.deseq<-function(data, design, designFormula, output="result") # output can also be "rlog" for rlog transformed counts
{
  data<-apply(data, c(1,2), round)
  colData<-data.frame(row.names = rownames(design))
  for (i in 1:dim(design)[2]) {
    colData[,i]<-as.factor(design[,i])
  }
  colnames(colData)<-colnames(design)
  if (missing(designFormula)){
    designFormula=paste0("~", colnames(colData))
  }
  dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design = formula(designFormula))
  if(dim(design)[2]==1 & length(unique(design[,1]))==2){
    dds <- DESeq(dds) # Wald test for pairwise
  } else {
    dds <- DESeq(dds, test="LRT", reduced=~1) # Likelihood ratio test for multiple groups/factors
  }
  if(output=="result"){
    res <- results(dds)
    res <- res[order(res$padj),]
    b<-match(rownames(res),rownames(data))
    out<-data.frame(geneID=rownames(res),order=b, res)
    out
  } else if(output=="rlog"){
    rld <- rlog(dds, blind=FALSE)
    data_tf<-as.data.frame(assay(rld))
    data_tf
  } else if(output=="vst"){
    vsd<-varianceStabilizingTransformation(dds, blind = FALSE)
    data_tf<-as.data.frame(assay(vsd))
    data_tf
  }
}

.bayseq<-function(data, replicates, design){
    groups <- list(NDE = c(rep(1,dim(design)[1])))
    groups[[colnames(design)[1]]]= design[,1]
    genenames<-rownames(data)
    data=as.matrix(data)
    rownames(data)<-1:length(rownames(data))
    cD <- new("countData", data = data, replicates = replicates, groups = groups)
    libsizes(cD) <- getLibsizes(cD)
    cD <- getPriors.NB(cD,  samplesize = 1000, cl = NULL)
    cD <- getLikelihoods(cD)
    cD@annotation <- data.frame(annotation = genenames)
    res<-topCounts(cD,group= colnames(design)[1],number=nrow(data))
    sizes<-as.numeric(libsizes(cD)/median(libsizes(cD)))
    baseMean<-rowMeans(t(t(data)/sizes))
    class<-unique(design[,1])
    log2FC<-rowMeans(t(t(data[,which(design[,1]==class[2])])/sizes[which(design[,1]==class[2])]))/rowMeans(t(t(data[,which(design[,1]==class[1])])/sizes[which(design[,1]==class[1])])) #reports log2FC of the first two classes by default
    out<-data.frame(geneID=res$annotation,order=rownames(res),baseMean=baseMean,log2FC=log2FC,padj=res[,dim(res)[2]-1], res[,c((dim(res)[2]-3):dim(res)[2])])
    out
}

.edger<-function(data, design, designFormula, output="result")
{
  colData<-data.frame(row.names = rownames(design))
  for (i in 1:dim(design)[2]) {
    colData[,i]<-as.factor(design[,i])
  }
  colnames(colData)<-colnames(design)
  if (missing(designFormula)) {
    designFormulaMatrix<-model.matrix(formula(paste0("~", colnames(colData))), colData)
  } else {
    designFormulaMatrix<-model.matrix(formula(designFormula), colData)
  }
  y <- DGEList(data)
  y <- calcNormFactors(y, method = "TMMwzp")
  if(output=="result"){
    if(dim(design)[2]==1 & length(unique(design[,1]))==2){ # exact test for pairwise
      y <- estimateDisp(y, designFormulaMatrix)
      fit <- glmQLFit(y,designFormulaMatrix)
      qlf <- glmQLFTest(fit,coef=2)
      qlf <- qlf$table
    } else {
      y <- estimateDisp(y, designFormulaMatrix) # GLM for multi-group/factor
      fit <- glmQLFit(y, designFormulaMatrix)
      qlf <- glmQLFTest(fit, coef=2:length(colnames(designFormulaMatrix)))
      qlf <- qlf$table
    }
    padj<-p.adjust(qlf$PValue,method='BH')
    qlf<-data.frame(qlf,padj)
    oe <- order(qlf[,"padj"])
    res<-qlf[oe,]
    baseMean<-rowMeans(t(t(data)/y$samples$norm.factors))
    b<-match(rownames(res),rownames(data))
    out<-data.frame(geneID=rownames(res),order=b,baseMean=baseMean[b], res)
    out
  } else if (output=="tmm"){
    tmm <- as.data.frame(cpm(y))
    tmm
  }
}

.save_pheatmap_pdf <- function(heatmap, filename) {
  pdf(filename)
  grid::grid.newpage()
  grid::grid.draw(heatmap$gtable)
  dev.off()
}

.FilterByPAdj <- function(genesDE,p=0.05, n=NULL)
{
  if(is.null(n)){
    ord <- which(genesDE$padj<p)
    as.character(genesDE$geneID[ord])
} else {
    ord <- order(genesDE$padj)
    as.character(genesDE$geneID[ord[1:n] ])
}
}

.ggplot2Col <- function(n){
  h = c(15,375)
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


.PlotCounts <- function(totalCounts, design){
  if(rownames(totalCounts)==rownames(design)){
    totalCounts$samples<-factor(rownames(totalCounts), rownames(totalCounts))
    design$samples<-factor(rownames(design), rownames(design))
    CountsPlots<-merge(totalCounts, design, by="samples")
  }
  CountsPlots_figures<-list()
  CountsPlots_figures[[1]]<-ggplot(CountsPlots, aes(y=genomeCounts, x=samples,fill=as.character(CountsPlots[,5])))+
    geom_bar(stat="summary", , funy='mean', width=0.3)+
    coord_cartesian(ylim = c(0, 1.2*max(CountsPlots$genomeCounts))) + theme_classic()+
    xlab("Samples")+ ylab("Mean genome-mapping spacer counts")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  CountsPlots_figures[[2]]<-ggplot(CountsPlots, aes(y=plasmidCounts, x=samples,fill=as.character(CountsPlots[,5])))+
    geom_bar(stat="summary", , funy='mean', width=0.3)+
    coord_cartesian(ylim = c(0, 1.2*max(CountsPlots$plasmidCounts))) + theme_classic()+
    xlab("Samples")+ ylab("Mean plasmid-mapping spacer counts")+theme(axis.text.x = element_text(angle = 90, hjust = 1))
  CountsPlots_figures
  ## For calculating mean counts and plotting with SEs ##
  # CountsPlots <-data.frame(unique(design))
  # colnames(CountsPlots)<-colnames(design)
  # MeanGenomeCounts<-c()
  # GenomeCountSEs<-c() #SE = standard error
  # MeanPlasmidCounts<-c()
  # PlasmidCountSEs<-c()
  # 
  # for(k in 1:dim(CountsPlots)[1]) {
  #   MeanGenomeCounts<-c(MeanGenomeCounts, mean(totalCounts$genomeCounts[which(design[,1]==CountsPlots[k,1])]))
  #   GenomeCountSEs<-c(GenomeCountSEs, sd(totalCounts$genomeCounts[which(design[,1]==CountsPlots[k,1])])/sqrt(length(which(design[,1]==CountsPlots[k,1]))))
  #   MeanPlasmidCounts<-c(MeanPlasmidCounts, mean(totalCounts$plasmidCounts[which(design[,1]==CountsPlots[k,1])]))
  #   PlasmidCountSEs<-c(PlasmidCountSEs, sd(totalCounts$plasmidCounts[which(design[,1]==CountsPlots[k,1])])/sqrt(length(which(design[,1]==CountsPlots[k,1]))))
  #   }
  # 
  # CountsPlots<-cbind(CountsPlots, MeanGenomeCounts, MeanPlasmidCounts, GenomeCountSEs, PlasmidCountSEs)
  # CountsPlots[is.na(CountsPlots)]=0
  # 
  # CountsPlots_figures<-list()
  # CountsPlots_figures[[1]]<-ggplot(CountsPlots, aes(y=MeanGenomeCounts, x=as.character(CountsPlots[,1])))+
  #   geom_bar(stat="identity", width=0.3, fill="deepskyblue3")+
  #   coord_cartesian(ylim = c(0, 1.5*max(CountsPlots$MeanGenomeCounts))) + theme_classic()+
  #   geom_errorbar(ymin=MeanGenomeCounts-GenomeCountSEs, ymax=MeanGenomeCounts+GenomeCountSEs,linetype=5, width = 0.1, color="darkblue")+
  #   xlab(colnames(design)[1])+ ylab("Mean genome-mapping spacer counts")
  # CountsPlots_figures[[2]]<-ggplot(CountsPlots, aes(y=MeanPlasmidCounts, x=as.character(CountsPlots[,1])))+
  #   geom_bar(stat="identity", width=0.3, fill="deepskyblue3")+
  #   coord_cartesian(ylim = c(0, 1.5*max(CountsPlots$MeanPlasmidCounts))) + theme_classic()+ ylab("Mean plasmid-mapping spacer counts")+
  #   geom_errorbar(ymin=MeanPlasmidCounts-PlasmidCountSEs, ymax=MeanPlasmidCounts+PlasmidCountSEs,linetype=5, width = 0.1, color="darkblue" )+
  #   xlab(colnames(design)[1])
  # CountsPlots_figures
}

.DEStats <- function(DEList, outPath, designFormula, K, vennDiagrams, geneBoxAndWhiskerPlots=TRUE, totalCountsFile, PCAplots=TRUE, clustering=TRUE, no=no, transformation) {         # formula to be used for Differential Expression.
  #If designFormula is null, all columns in the design matrix will be individually used as classes.
  #If supplied, only this formula will be used

  data<-as.data.frame(DEList[[1]])
  design<-as.data.frame(DEList[[2]])
  if(!is.null(totalCountsFile)){
    totalCounts<-as.data.frame(DEList[[3]])
  }
  

  union_de <- list()
  union_de_pval <-list()
  intersect_DE <- list()

  # creating a vector of replicates

  replicates<-c()
  rep=1
  for(rw in 1:dim(design)[1]){
    if(rw>1) {
      if(isTRUE(all.equal(design[rw,],design[rw-1,],check.attributes = FALSE)))
      {
        replicates<-c(replicates,rep)
      } else {
        rep=rep+1
        replicates<-c(replicates,rep)
      }
    } else {
      replicates<-c(replicates,rep)
    }
  }

  ## SCRIPT PARSES THROUGH EACH COLUMN OF DESIGN MATRIX ##

  for(i in 1:length(colnames(design))) {

    if (is.null(transformation)){
      data_transformed <- data
    } else if (transformation=="log2"){
      # log2 table of counts data
      data_transformed <- as.data.frame(apply(data, c(1,2), function(x) log2(x+1)))
    } else if (transformation=="rlog"){
      data_transformed <- .deseq(data, design[,i, drop=FALSE], output = "rlog")
    } else if (transformation=="vst"){
      data_transformed <- .deseq(data, design[,i, drop=FALSE], output = "vst")
    } else if (transformation=="tmm"){
      data_transformed <- .edger(data, design[,i, drop=FALSE], output = "tmm")
      data_transformed <- as.data.frame(apply(data_transformed, c(1,2), function(x) log2(x+1)))
    } else if (transformation=="thinCounts"){
      minLibSize <- min(colSums(data, na.rm = FALSE, dims = 1), na.rm = FALSE)
      data_transformed <- thinCounts(data, prob=NULL, target.size=minLibSize)
      data_transformed <- as.data.frame(apply(data_transformed, c(1,2), function(x) log2(x+1)))
    }


    ## PLOTTING TOTAL COUNTS BY SAMPLE ##
    if(!is.null(totalCountsFile)){
      countplots <- .PlotCounts(totalCounts, design[,i, drop=FALSE])
      pdf(paste0(outPath, "/totalCountsBySamplePlots_",colnames(design)[i], ".pdf"))
      for(l in 1:2){print(countplots[[l]])}
      dev.off()
    }

    ## Box and Whisker plots of transformed Spacer counts by gene for each sample ##
    if(geneBoxAndWhiskerPlots) {
      geneSort<-order(-rowSums(data_transformed))
      GeneBoxPlots<-t(data_transformed[geneSort[1:no],])
      GeneBoxPlots<-melt(GeneBoxPlots)
      colnames(GeneBoxPlots)<-c("SampleID", "Gene", "transformed_SpacerCounts")
      GeneBoxPlots[,4]<-0
      colnames(GeneBoxPlots)[4]<-"design"
      for(t in 1:dim(data)[2]) {
        GeneBoxPlots[t,4]<-as.character(design[which(rownames(design)==GeneBoxPlots$SampleID[t]), i])
      }
      GeneBoxPlots[(dim(data)[2]+1):dim(GeneBoxPlots)[1],4]<-as.character(rep(GeneBoxPlots[1:dim(data)[2],4], (dim(GeneBoxPlots)[1]/dim(data)[2])-1))
      GeneBoxPlots$SampleIDs<-as.character(rep(paste0('S',1:dim(data)[2]), dim(GeneBoxPlots)[1]/dim(data)[2]))
      GeneBoxPlots$SampleIDs<-factor(GeneBoxPlots$SampleID, GeneBoxPlots$SampleID[1:dim(data)[2]])
      ggplot(data=GeneBoxPlots, aes(x=SampleIDs,y=transformed_SpacerCounts))+
        geom_boxplot(aes(fill=design))+theme_classic()+
        xlab("Sample ID") + ylab(paste0(as.character(transformation), " transformed gene-mapping spacer counts"))+ggtitle(paste0(as.character(transformation), " transformed gene-mapping spacer counts for top ", as.character(no), " genes"))+theme(axis.text.x = element_text(angle = 90, hjust = 1))
      ggsave(paste0(outPath, "/GeneBoxPlots_", colnames(design)[i],".pdf"), height = 8.5, width = 10)
    }

    # PCA ##
    if(PCAplots) {
      if(dim(data)[1]<no) {
        no=dim(data)[1]
      }
      sds <- apply(data_transformed, 1, sd)
      o <- order(sds, decreasing = TRUE)
      all_pca<-as.data.frame(t(data_transformed[o[1:no],]))
      all_pca[,no+1]<-as.character(design[,i])
      colnames(all_pca)[no+1]<-colnames(design)[i]
      colnames(all_pca)[1:no]<-1:no
      autoplot(prcomp(all_pca[,1:no]), data = all_pca, size=4, colour = colnames(all_pca)[no+1]) + theme_bw()+ ggtitle(paste0("PCA plot for top ", as.character(no)," genes sorted by variance")) + scale_size(guide="none")
      ggsave(paste0(outPath, "/PCA_", colnames(design)[i],".svg"), height = 8.5, width = 10)
      if(clustering){
        autoplot(fanny(all_pca[,1:no], K), data = all_pca, size=4,shape = colnames(all_pca)[no+1], frame = TRUE, frame.type = 'norm') + theme_bw() + ggtitle(paste0("FANNY clustering for top ", as.character(no)," genes sorted by variance")) + scale_size(guide="none")
        ggsave(paste0(outPath, "/FANNY_", colnames(design)[i],".pdf"), height = 8.5, width = 10)
        set.seed(1)
        autoplot(kmeans(all_pca[,1:no], K), data = all_pca, size=4, shape = colnames(all_pca)[no+1], frame = TRUE, frame.type = 'norm') +  theme_bw() + ggtitle(paste0("K-means clustering for top ", as.character(no)," genes sorted by variance")) + scale_size(guide="none")
        ggsave(paste0(outPath, "/KMeans_", colnames(design)[i],".pdf"), height = 8.5, width = 10)
      }
    }

    ## DIFFERENTIAL EXPRESSION ##
    if (is.null(designFormula)|is.na(designFormula)) {
      out.de <- .deseq(data, design[,i, drop=FALSE])
      out.er <- .edger(data, design[,i, drop=FALSE])
      out.bs <- .bayseq(data, replicates, design[,i, drop=FALSE])

      # rownames in the result tables as IDs - for consistency
      rownames(out.de) <- out.de$geneID
      rownames(out.er) <- out.er$geneID
      rownames(out.bs) <- out.bs$geneID

      # finding top genes
      top.de <- .FilterByPAdj(out.de, n=20)
      top.er <- .FilterByPAdj(out.er, n=20)
      top.bs <- .FilterByPAdj(out.bs, n=20)
      top.de_pval <- .FilterByPAdj(out.de)
      top.er_pval <- .FilterByPAdj(out.er)
      top.bs_pval <- .FilterByPAdj(out.bs)

      # # union of top 20 genes from DE tools
      union_de[[colnames(design)[i]]] <- union(top.de, union(top.er, top.bs))
      # union of top genes by p value from DE tools
      union_de_pval[[colnames(design)[i]]] <- union(top.de_pval, union(top.er_pval, top.bs_pval))
      # intersect of top genes from DE tools
      intersect_DE[[colnames(design)[i]]] <- intersect(top.de_pval, intersect(top.er_pval, top.bs_pval))

      ## Venn Diagram ##
      if (vennDiagrams){

        pdf(paste0(outPath, "/", "VENN_", colnames(design)[i],".pdf"), width=10, height=10)
        draw.triple.venn(area1 = length(top.de_pval), area2 = length(top.er_pval), area3 = length(top.bs_pval), n12 = length(intersect(top.de_pval,top.er_pval)), n23 = length(intersect(top.er_pval,top.bs_pval)), n13 = length(intersect(top.bs_pval,top.de_pval)),
                         n123 = length(intersect(top.de_pval,intersect(top.er_pval,top.bs_pval))), category = c("DESeq2", "edgeR", "baySeq"),
                         col = FALSE, fill = c("deepskyblue2", "red", "green"), alpha = 0.3)
        dev.off()
      }


      ## Heatmaps ##

      annot_col<-data.frame(as.character(design[,i]))
      rownames(annot_col)=rownames(design)
      colnames(annot_col)=colnames(design)[i]
      colors<-c()
      for(j in 1:length(design[,i])){
        colors<-c(colors, .ggplot2Col(length(unique(design[,i])))[which(unique(design[,i])==design[j,i])])
      }
      colors_corr <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
      out_union <-  pheatmap(data_transformed[union_de_pval[[colnames(design)[i]]], ], annotation_col = annot_col, fontsize = 8, fontsize_row = 5, fontsize_col = 7, main = "heatmap for union of DE genes from all tools")
      if(length(intersect_DE[[colnames(design)[i]]])>2) {
        out_intersect<-pheatmap(data_transformed[intersect_DE[[colnames(design)[i]]], ], annotation_col = annot_col, fontsize = 8, fontsize_row = 5, fontsize_col = 7, main = "heatmap for intersect of DE genes from all tools")
      }

      if(length(intersect_DE[[colnames(design)[i]]])>2){
        heatmap<-pheatmap(data_transformed[intersect_DE[[colnames(design)[i]]], ], annotation_col = annot_col, fontsize = 8, fontsize_row = 5, fontsize_col = 7, main = "intersect of significant DE genes (adjusted p-val < 0.05) from all tools")
        .save_pheatmap_pdf(heatmap, paste0(outPath, "/", "heatmaps_intersect_", colnames(design)[i],".pdf"))
        }
      heatmap<-pheatmap(data_transformed[union_de[[colnames(design)[i]]], ], annotation_col = annot_col, fontsize = 8, fontsize_row = 5, fontsize_col = 7, main = "union of top 20 DE genes from all tools")
      .save_pheatmap_pdf(heatmap, paste0(outPath, "/", "heatmaps_union_top20_", colnames(design)[i],".pdf"))
      heatmap<-pheatmap(data_transformed[union_de_pval[[colnames(design)[i]]], ], annotation_col = annot_col, fontsize = 8, fontsize_row = 5, fontsize_col = 7, main = "union of significant DE genes (adjusted p-val < 0.05) from all tools")
      .save_pheatmap_pdf(heatmap, paste0(outPath, "/", "heatmaps_union_", colnames(design)[i],".pdf"))
      heatmap<-pheatmap(data_transformed[o[1:50], ], annotation_col = annot_col, fontsize = 8, fontsize_row = 5, fontsize_col = 7, main = "top 50 genes sorted by variance")
      .save_pheatmap_pdf(heatmap, paste0(outPath, "/", "heatmaps_topvariance_", colnames(design)[i],".pdf"))
      pdf(paste0(outPath, "/", "heatmaps_correlation_", colnames(design)[i],".pdf"))
        heatmap.2(cor(data), trace="none",RowSideColors = as.character(colors), col=colors_corr, ColSideColors = as.character(colors), margins = c(10,10))
        legend("topright", inset=.02, title=colnames(design[i]), legend=unique(design[,i]), fill=.ggplot2Col(length(unique(design[,i]))), cex=0.5)
      dev.off()



      ## write the differential expression results to an output file ##
      out_union<-data.frame(geneID=out_union$tree_row[["labels"]][out_union$tree_row[["order"]]], padj_DESeq2=NA, padj_edgeR=NA, padj_bayseq=NA)
      for(ind in 1:dim(out_union)[1]){
        if(out_union$geneID[ind]%in%out.de$geneID){
          out_union$padj_DESeq2[ind]=out.de$padj[which(out.de$geneID==as.character(out_union$geneID[ind]))]
        }
        if(out_union$geneID[ind]%in%out.er$geneID){
          out_union$padj_edgeR[ind]=out.er$padj[which(out.er$geneID==as.character(out_union$geneID[ind]))]
        }
        if(out_union$geneID[ind]%in%out.bs$geneID){
          out_union$padj_bayseq[ind]=out.bs$padj[which(out.bs$geneID==as.character(out_union$geneID[ind]))]
        }
      }
      if(length(intersect_DE[[colnames(design)[i]]])>2){
        out_intersect<-data.frame(geneID=out_intersect$tree_row[["labels"]][out_intersect$tree_row[["order"]]], padj_DESeq2=NA, padj_edgeR=NA, padj_bayseq=NA)
        for(ind in 1:dim(out_intersect)[1]){
          out_intersect$padj_DESeq2[ind]=out.de$padj[which(out.de$geneID==as.character(out_intersect$geneID[ind]))]
          out_intersect$padj_edgeR[ind]=out.er$padj[which(out.er$geneID==as.character(out_intersect$geneID[ind]))]
          out_intersect$padj_bayseq[ind]=out.bs$padj[which(out.bs$geneID==as.character(out_intersect$geneID[ind]))]
        }
        write.csv(out_intersect,   file=paste0(outPath, "/DEgenes_",colnames(design)[i],"_Intersect_SignatureGenes.csv"))
      }

      write.csv(out.de, file=paste0(outPath, "/", "DEseq_", colnames(design)[i],".csv"))
      write.csv(out.bs, file=paste0(outPath, "/", "baySeq_", colnames(design)[i],".csv"))
      write.csv(out.er, file=paste0(outPath, "/", "edgeR_", colnames(design)[i],".csv"))
      write.csv(out_union,   file=paste0(outPath, "/DEgenes_",colnames(design)[i],"_Union.csv"))
    }
  }

  if (!is.null(designFormula)&!is.na(designFormula)){
    if(!startsWith(designFormula,"~")) {
      stop("Design formula must start with '~' symbol!")
    }
    if(!isEmpty(grep(":",designFormula))){
      Ratios<-list()
    }
    tempDesign<-strsplit2(designFormula,"+", fixed=TRUE)
    tempDesign[1,1]<-gsub("~", "", tempDesign[1,1], fixed = TRUE)
    SingleElements<-c()
    for(i in 1:dim(tempDesign)[2]) {
      if(isEmpty(grep(":",tempDesign[1,i]))) {
        if(isEmpty(grep(tempDesign[1,i], colnames(design)))) {
          stop("Check design Formula! Only use EXACT colnames from design matrix and don't use symbols except + and :")
        } else {
          SingleElements<-c(SingleElements,tempDesign[1,i])
        }
      } else {
        tempRatioElements<-strsplit2(tempDesign[1,i],":", fixed=TRUE)
        for(m in 1:dim(tempRatioElements)[2]) {
          if(isEmpty(grep(tempRatioElements[1,m], colnames(design)))) {
            stop("Check design Formula! Only use EXACT colnames from design matrix and don't use symbols except + and :")
          } else {
            SingleElements<-c(SingleElements,tempRatioElements[1,m])
          }
        }
      }
      SingleElements<-unique(SingleElements)
    }
    out.de <- .deseq(data, design[,which(colnames(design)%in%SingleElements), drop=FALSE], designFormula)
    out.er <- .edger(data, design[,which(colnames(design)%in%SingleElements), drop=FALSE], designFormula)

    # finding top genes
    top.de <- .FilterByPAdj(out.de, n=30)
    top.er <- .FilterByPAdj(out.er, n=30)
    top.de_pval <- .FilterByPAdj(out.de)
    top.er_pval <- .FilterByPAdj(out.er)

    # union of top genes from DE tools
    union_de <- union(top.de, top.er)
    # intersect of top genes from DE tools
    intersect_DE <- intersect(top.de_pval, top.er_pval)
    union_de_pval <- union(top.er_pval, top.de_pval)

    ## Venn Diagram ##
    if (vennDiagrams){
      pdf(paste0(outPath, "/", "VENN_", designFormula,".pdf"), width=10, height=10)
      draw.pairwise.venn(area1 = length(top.de_pval), area2 = length(top.er_pval),  cross.area = length(intersect(top.de_pval,top.er_pval)), category = c("DESeq2", "edgeR"),
                         col = c("deepskyblue2", "red"))
      dev.off()
    }

    ## Heatmaps ##

    ## Heatmaps ##

    annot_col<-as.data.frame(design[,which(colnames(design)%in%SingleElements), drop=FALSE])
    colors_corr <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
    out_union <-  pheatmap(data_transformed[union_de, ], annotation_col = annot_col, fontsize = 8, fontsize_row = 5, fontsize_col = 7, main = "heatmap for union of DE genes from all tools")
    if(length(intersect_DE)>2) {
      out_intersect<-pheatmap(data_transformed[intersect_DE, ], annotation_col = annot_col, fontsize = 8, fontsize_row = 5, fontsize_col = 7, main = "heatmap for intersect of DE genes from all tools")
    }

    if(length(intersect_DE)>2){
      heatmap<-pheatmap(data_transformed[intersect_DE, ], annotation_col = annot_col, fontsize = 8, fontsize_row = 5, fontsize_col = 7, main = "intersect of significant DE genes (adjusted p-val < 0.05) from all tools")
      .save_pheatmap_pdf(heatmap, paste0(outPath, "/", "heatmaps_intersect_", designFormula,".pdf"))
    }
    heatmap<-pheatmap(data_transformed[union_de, ], annotation_col = annot_col, fontsize = 8, fontsize_row = 5, fontsize_col = 7, main = "union of top 20 DE genes from all tools")
    .save_pheatmap_pdf(heatmap, paste0(outPath, "/", "heatmaps_union_top20_", designFormula,".pdf"))
    heatmap<-pheatmap(data_transformed[union_de_pval, ], annotation_col = annot_col, fontsize = 8, fontsize_row = 5, fontsize_col = 7, main = "union of significant DE genes (adjusted p-val < 0.05) from all tools")
    .save_pheatmap_pdf(heatmap, paste0(outPath, "/", "heatmaps_union_", designFormula,".pdf"))
    heatmap<-pheatmap(data_transformed[o[1:50], ], annotation_col = annot_col, fontsize = 8, fontsize_row = 5, fontsize_col = 7, main = "top 50 genes sorted by variance")
    .save_pheatmap_pdf(heatmap, paste0(outPath, "/", "heatmaps_topvariance_", designFormula,".pdf"))
    pdf(paste0(outPath, "/", "heatmaps_correlation_", designFormula,".pdf"))
    heatmap.2(cor(data), trace="none", col=colors_corr,
              margins = c(10,10))
    dev.off()

    ## write the differential expression results to an output file ##

    out_union<-data.frame(geneID=out_union$tree_row[["labels"]][out_union$tree_row[["order"]]], padj_DESeq2=NA, padj_edgeR=NA)
    for(ind in 1:dim(out_union)[1]){
      if(out_union$geneID[ind]%in%out.de$geneID){
        out_union$padj_DESeq2[ind]=out.de$padj[which(out_union$geneID[ind]==out.de$geneID)]
      }
      if(out_union$geneID[ind]%in%out.er$geneID){
        out_union$padj_edgeR[ind]=out.er$padj[which(out_union$geneID[ind]==out.er$geneID)]
      }
    }
    if(length(intersect_DE)>2){
      out_intersect<-data.frame(geneID=out_intersect$tree_row[["labels"]][out_intersect$tree_row[["order"]]], padj_DESeq2=NA, padj_edgeR=NA)
      for(ind in 1:dim(out_intersect)[1]){
        out_intersect$padj_DESeq2[ind]=out.de$padj[which(out_intersect$geneID[ind]==out.de$geneID)]
        out_intersect$padj_edgeR[ind]=out.er$padj[which(out_intersect$geneID[ind]==out.er$geneID)]
      }
      write.csv(out_intersect,   file=paste0(outPath, "/DEgenes_",designFormula,"_IntersectSignature.csv"))
    }

    write.csv(out.de, file=paste0(outPath, "/", "DEseq_", designFormula,".csv"))
    write.csv(out.er, file=paste0(outPath, "/", "edgeR_", designFormula,".csv"))
    write.csv(out_union,   file=paste0(outPath, "/DEgenes_",designFormula,"_Union.csv"))
  }
}



