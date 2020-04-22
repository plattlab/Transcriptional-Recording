# Transcriptional Recordings using Record-seq

This repository contains python scripts and a R package for analysing datasets generated using the Record-seq workflow described at https://www.nature.com/articles/s41586-018-0569-1

## Getting Started

Download the scripts within the ```primary-analysis``` directory and save the ```spacerExtractor.py``` and ```SErmdup.py``` scripts in your preferred scripts directory. Create a dedicated directory for the analysis of your experiment at your preferred path - by either typing ```$ mkdir experimentDirectory``` in a bash session or using your UI. Download and save the ```Snakefile``` and ```config.yml``` file within this directory. 

Further, from within a session of R \(using RStudio or an R session within bash \), download the ```recoRdseq``` R package using:

```
> install.packages("devtools")
> library(devtools)
> install_github("plattlab/Transcriptional-Recording", subdir="recoRdseq")
```

## License

This repository is licensed under the Apache License Version 2.0

## Authors

* **Tanmay Tanna** - (https://github.com/TanmayTanna)

