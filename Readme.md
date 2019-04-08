# Transcriptional Recordings using Record-seq

This repository contains python scripts and a R package for analysing datasets generated using the Record-seq workflow described at https://www.nature.com/articles/s41586-018-0569-1

## Getting Started

Download the scripts within the ```primary-analysis``` directory and save the ```spacerExtractor.py``` and ```SErmdup.py``` scripts in your preferred directory. Create a dedicated directory for the analysis of your experiment either typing ```$ mkdir experimentDirectory``` in a bash session or using your UI at your preferred path and save the ```Snakefile``` and ```config.yml``` file within this directory. 

Further, from within a session of R \(using RStudio or an R session within bash \), download the ```recoRdseq``` R package using:

```
> install.packages("devtools")
> library(devtools)
> install_github("plattlab/Transcriptional-Recording", subdir="recoRdseq")
```



### Installing prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
