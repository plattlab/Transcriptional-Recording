\name{recordSeqAnalysis}
\alias{recordSeqAnalysis}

\title{ recordSeqAnalysis
}
\description{
Does all the analyses and plots for recoRdseq
#\code{seq\_read} table
}
\usage{
recordSeqAnalysis (
  FeatureCounts = "../data_Tanmay_paraquat/UMIGenomeCounts2.txt",
  designMatrix = "../data_Tanmay_paraquat/DesignMatrix_timepoint3.xlsx",
  TotalCountsbySample = "../data_Tanmay_paraquat/TotalCountsbySample.txt",
  outPath = "./output",
  designFormula =  NA,
  mil=1000000,
  MinCountsPerSample = 100,
  GeneBoxAndWhiskerPlots = TRUE,
  KMeans_Cluster = TRUE,
  CountBarPlots = TRUE,
  PCAplots = TRUE,
  VennDiagrams =TRUE,
  ProcessInput = TRUE)
}
\arguments{
  \item{FeatureCounts}{Input count table
}
  \item{designMatrix}{Experimental design matrix
}
  \item{TotalCountsbySample}{TotalCountsbySample
}
  \item{outPath}{Output directory
}
  \item{designFormula}{designFormula
}
  \item{mil}{mil
}
  \item{MinCountsPerSample}{Threshold for counts
}
  \item{GeneBoxAndWhiskerPlots}{GeneBoxAndWhiskerPlots
}
  \item{KMeans_Cluster}{KMeans_Cluster
}
  \item{CountBarPlots}{CountBarPlots
}
  \item{PCAplots}{PCAplots
}
  \item{VennDiagrams}{VennDiagrams
}
  \item{ProcessInput}{ProcessInput
}

\value{ None
}
\author{
   Michal Okoniewski, Tanmay Tanna
}
\examples{
#recordSeqAnalysis (
#  FeatureCounts = "../data_Tanmay_paraquat/UMIGenomeCounts2.txt",    #count table
#  designMatrix = "../data_Tanmay_paraquat/DesignMatrix_timepoint3.xlsx",
# TotalCountsbySample = "../data_Tanmay_paraquat/TotalCountsbySample.txt",
# outPath = "./output")
}
