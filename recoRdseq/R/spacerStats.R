library(ggplot2)
library(readr)
library(readxl)

spacerStats<- function(
  summaryStats = NULL,
  SpacerInfoFiles_Path=NULL, # directory containing spacer info files generated during spacer extraction
  outPath="/path/for/output",
  designMatrix = "path/to/designMatrix"
  )
{
  if(!dir.exists(outPath)){dir.create(outPath)}

  design <- as.data.frame(read_excel(designMatrix))
  rownames(design) <- design[,1]
  design <- design[,-1, drop=FALSE]

  ## Reading in stats files and pre-processing them

  if(!is.null(summaryStats)){
    stats<-as.data.frame(read.delim(summaryStats, header = TRUE, sep="\t"))
    rownames(stats)<-stats[,1]
    stats<-stats[,-1]
    stats<-stats[which(rownames(stats)%in%rownames(design)),]
    stats<-stats[match(rownames(design),rownames(stats)),]
    stats$uniqueSpacersPerMillionReads = stats$uniqueSpacers*1000000/stats$totalReads

    ## Looping through the columns of the design matrix and creating spacer plots using summaryStats
    for(i in 1:length(colnames(design))) {
        spacerplots <- .SpacerPlot(stats, design[,i, drop=FALSE])
        pdf(paste0(outPath, "/summarySpacerPlots_",colnames(design)[i], ".pdf"))
        for(l in 1:4){print(spacerplots[[l]])}
        dev.off()
    }
  }
  ## Creating strings for info file paths and reading in info files
  if(!is.null(SpacerInfoFiles_Path)){
    summary_filepaths=list.files(path = SpacerInfoFiles_Path, pattern = "*.info.txt")
    summary_files<-list()
    
    for(i in 1: length(summary_filepaths)){
      sampleName=gsub(".info.txt", "", summary_filepaths[i])
      summary_files[[sampleName]]<- read_delim(paste0(SpacerInfoFiles_Path, "/", summary_filepaths[i]),
                                               "\t", escape_double = FALSE, trim_ws = TRUE)
    }
    
    ## Making GC content and spacer length plots for spacer files ##
    
    for(i in 1:length(colnames(design))) {
      
      
      designFactors<-unique(design[,i])
      
      for(df in 1:length(designFactors)){
        GC_content<-c()
        Sequence_Length<-c()
        samples<-rownames(design)[which(design[,i]==designFactors[df])]
        
        for(s in 1:length(samples)){
          GC_content<-c(GC_content,summary_files[[samples[s]]]$GC_content)
          Sequence_Length<-c(Sequence_Length, summary_files[[samples[s]]]$Sequence_Length)
        }
        
        GC_content<-as.data.frame(GC_content)
        Sequence_Length<-as.data.frame(Sequence_Length)
        Sequence_Length_freq<-data.frame(table(Sequence_Length))
        k<-sum(Sequence_Length_freq$Freq)
        Sequence_Length_freq$Freq<-Sequence_Length_freq$Freq*100/k
        
        GC_content_distribution<-data.frame(GC_content=seq(0,99,by=5), frequency=0)
        for(j in 1:20){GC_content_distribution$frequency[j]=as.numeric(length(which(GC_content$GC_content>GC_content_distribution$GC_content[j]&GC_content$GC_content<(GC_content_distribution$GC_content[j]+5))))}
        k = sum(GC_content_distribution$frequency)
        GC_content_distribution$frequency<-GC_content_distribution$frequency*100/k
        
        ggplot(GC_content_distribution, aes(x=GC_content, y=frequency))+
          geom_bar(stat="identity", fill="deepskyblue3")+
          theme_classic(base_size = 7)+
          scale_x_continuous(breaks = seq(0, 100, by = 10))+
          xlab("GC content (%)")+
          ylab("Frequency (%)")+
          geom_vline(aes(xintercept=50), linetype=3)
        ggsave(paste0(outPath,"/SpacerGCcontentDistribution_",colnames(design)[i],"_", designFactors[df], ".pdf"),  height = 8.5, width = 10)
        
        
        ggplot(Sequence_Length_freq, aes(x=Sequence_Length, y=Freq))+
          geom_bar(stat="identity", fill="deepskyblue3")+
          theme_classic(base_size = 7)+
          xlab("Spacer Length")+
          ylab("Frequency (%)")
        ggsave(paste0(outPath,"/SpacerLengthFrequencies_",colnames(design)[i],"_", designFactors[df], ".pdf"),  height = 8.5, width = 10)
      }
    }
  }
}


  ## Function to create spacer plots

.SpacerPlot <- function(stats, design){
  SpacerPlots <-data.frame(unique(design))
  colnames(SpacerPlots)<-colnames(design)
  MeanUniqueSpacers<-c()
  MeanUniqueSpacersPerMillionReads<-c()
  UniqueSpacerSEs<-c()
  UniqueSpacersPerMillionReads_SEs<-c()
  MeanUniqueSingleAcquisitions<-c()
  UniqueSingleAcquisitionSEs<-c()
  MeanUniqueDoubleAcquisitions<-c()
  UniqueDoubleAcquisitionSEs<-c()
  for(k in 1:dim(SpacerPlots)[1]) {
    MeanUniqueSpacers<-c(MeanUniqueSpacers, mean(stats$uniqueSpacers[which(design[,1]==SpacerPlots[k,1])]))
    MeanUniqueSpacersPerMillionReads<-c(MeanUniqueSpacersPerMillionReads, mean(stats$uniqueSpacersPerMillionReads[which(design[,1]==SpacerPlots[k,1])]))
    UniqueSpacerSEs<-c(UniqueSpacerSEs, sd(stats$uniqueSpacers[which(design[,1]==SpacerPlots[k,1])])/sqrt(dim(SpacerPlots)[1]))
    UniqueSpacersPerMillionReads_SEs<-c(UniqueSpacersPerMillionReads_SEs, sd(stats$uniqueSpacersPerMillionReads[which(design[,1]==SpacerPlots[k,1])])/sqrt(dim(SpacerPlots)[1]))
    MeanUniqueSingleAcquisitions<-c(MeanUniqueSingleAcquisitions, mean(stats$uniqueSingleAcquisitions[which(design[,1]==SpacerPlots[k,1])]))
    UniqueSingleAcquisitionSEs<-c(UniqueSingleAcquisitionSEs, sd(stats$uniqueSingleAcquisitions[which(design[,1]==SpacerPlots[k,1])])/sqrt(dim(SpacerPlots)[1]))
    MeanUniqueDoubleAcquisitions<-c(MeanUniqueDoubleAcquisitions, mean(stats$uniqueDoubleAcquisitions[which(design[,1]==SpacerPlots[k,1])]))
    UniqueDoubleAcquisitionSEs<-c(UniqueDoubleAcquisitionSEs, sd(stats$uniqueDoubleAcquisitions[which(design[,1]==SpacerPlots[k,1])])/sqrt(dim(SpacerPlots)[1]))
  }
  SpacerPlots$MeanUniqueSpacers<-MeanUniqueSpacers
  SpacerPlots$MeanUniqueSpacersPerMillionReads<-MeanUniqueSpacersPerMillionReads
  SpacerPlots$UniqueSpacerSEs<-UniqueSpacerSEs
  SpacerPlots$UniqueSpacersPerMillionReads_SEs<-UniqueSpacersPerMillionReads_SEs
  SpacerPlots$MeanUniqueSingleAcquisitions<-MeanUniqueSingleAcquisitions
  SpacerPlots$UniqueSingleAcquisitionSEs<-UniqueSingleAcquisitionSEs
  SpacerPlots$MeanUniqueDoubleAcquisitions<-MeanUniqueDoubleAcquisitions
  SpacerPlots$UniqueDoubleAcquisitionSEs<-UniqueDoubleAcquisitionSEs
  SpacerPlots[is.na(SpacerPlots)]=0
  

  spacerplots<-list()
  spacerplots[[1]]<-ggplot(SpacerPlots, aes(y=MeanUniqueSpacers, x=as.character(SpacerPlots[,1])))+
    geom_bar(stat="identity", width=0.3,  fill="deepskyblue3")+
    coord_cartesian(ylim = c(0, 1.5*max(SpacerPlots$MeanUniqueSpacers))) + theme_classic(base_size = 7)+
    geom_errorbar(ymin=MeanUniqueSpacers-UniqueSpacerSEs, ymax=MeanUniqueSpacers+UniqueSpacerSEs,linetype=5, width = 0.1, color="darkblue" )+
    xlab(colnames(design)[1]) + ylab("Mean Unique spacers")
  spacerplots[[2]]<-ggplot(SpacerPlots, aes(y=MeanUniqueSingleAcquisitions, x=as.character(SpacerPlots[,1])))+
    geom_bar(stat="identity", width=0.3,  fill="deepskyblue3")+
    coord_cartesian(ylim = c(0, 1.5*max(SpacerPlots$MeanUniqueSingleAcquisitions))) + theme_classic(base_size = 7)+
    geom_errorbar(ymin=MeanUniqueSingleAcquisitions-UniqueSingleAcquisitionSEs, ymax=MeanUniqueSingleAcquisitions+UniqueSingleAcquisitionSEs,linetype=5, width = 0.1, color="darkblue" )+
    xlab(colnames(design)[1]) + ylab("Mean unique single acquisitions")
  spacerplots[[3]]<-ggplot(SpacerPlots, aes(y=MeanUniqueDoubleAcquisitions, x=as.character(SpacerPlots[,1])))+
    geom_bar(stat="identity", width=0.3,  fill="deepskyblue3")+
    coord_cartesian(ylim = c(0, 1.5*max(SpacerPlots$MeanUniqueDoubleAcquisitions))) + theme_classic(base_size = 7)+
    geom_errorbar(ymin=MeanUniqueDoubleAcquisitions-UniqueDoubleAcquisitionSEs, ymax=MeanUniqueDoubleAcquisitions+UniqueDoubleAcquisitionSEs,linetype=5, width = 0.1, color="darkblue" )+
    xlab(colnames(design)[1]) + ylab("Mean Unique Double acquisitions")
  spacerplots[[4]]<-ggplot(SpacerPlots, aes(y=MeanUniqueSpacersPerMillionReads, x=as.character(SpacerPlots[,1])))+
    geom_bar(stat="identity", width=0.3,  fill="deepskyblue3")+
    coord_cartesian(ylim = c(0, 1.5*max(SpacerPlots$MeanUniqueSpacersPerMillionReads))) + theme_classic(base_size = 7)+
    geom_errorbar(ymin=MeanUniqueSpacersPerMillionReads-UniqueSpacersPerMillionReads_SEs, ymax=MeanUniqueSpacersPerMillionReads+UniqueSpacersPerMillionReads_SEs,linetype=5, width = 0.1, color="darkblue" )+
    xlab(colnames(design)[1]) + ylab("Mean Unique spacers per million sequencing reads")
  spacerplots
}

