library(ggplot2)
library(readr)
library(readxl)
input_path="./spacerStats" # directory containing spacer info files generated during spacer extraction
output_path="./spacerStats/plots"
if(!dir.exists(output_path)){dir.create(path = output_path)}


## Creating strings for info file paths and reading in info files 

length_filepaths=list.files(path = input_path, pattern = "*.spacerLengths.txt")
gc_filepaths=list.files(path = input_path, pattern = "*.spacerGC.txt")
gc_files<-list()

for(i in 1: length(gc_filepaths)){
  name=sub(".txt", "", gc_filepaths[i])
  gc<- as.data.frame(read_delim(paste0(input_path, "/", gc_filepaths[i]),
                                  "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
  
  gcFreq<-data.frame(GC_content=seq(0,99,by=5), frequency=0)
  for(i in 1:20){gcFreq$frequency[i]=sum(gc[which(gc[,1]>gcFreq$GC_content[i]&gc[,1]<gcFreq$GC_content[i+1]),2])}
  k = sum(gcFreq$frequency)
  gcFreq$frequency<-gcFreq$frequency*100/k
  
  ggplot(gcFreq, aes(x=GC_content, y=frequency))+
    geom_bar(stat="identity", fill="deepskyblue3")+
    theme_classic(base_size = 7)+
    scale_x_continuous(breaks = seq(0, 100, by = 10))+
    xlab("GC content (%)")+
    ylab("Frequency (%)")+
    geom_vline(aes(xintercept=50), linetype=3)
  ggsave(paste0(output_path,"/",name,".pdf"),  height = 8.5, width = 10)
}
  

length_files<-list()
for(i in 1: length(length_filepaths)){
  name=sub(".txt", "", length_filepaths[i])
  length<- as.data.frame(read_delim(paste0(input_path, "/", length_filepaths[i]),
                             "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE))
  
  colnames(length)<-c("Sequence_Length", "Freq")
  ggplot(length , aes(x=Sequence_Length, y=Freq))+
    geom_bar(stat="identity", fill="deepskyblue3")+
    theme_classic(base_size = 7)+
    xlab("Spacer Length")+
    ylab("Frequency (%)")
  ggsave(paste0(output_path,"/",name,".pdf"),  height = 8.5, width = 10)
}
