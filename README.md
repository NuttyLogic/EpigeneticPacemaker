# Epigenetic Pacemaker
## A fast conditional expectation maximization algorithm for modeling epigenetic state

DNA methylation is widely used to model physiological phenotypes, such as 
aging[1](https://doi.org/10.1186/gb-2013-14-10-r115) and type II diabetes[2](https://doi.org/10.1093/hmg/ddy093). 
The epigenetic pacemaker, **EPM**, is an implementation of a fast conditional expectation maximization algorithm used to 
model epigenetic states associated with a phenotype of interest [3](https://doi.org/10.2217/epi-2017-0130) The EPM was first introduced by Snir et al. 
[4](https://doi.org/10.1371/journal.pcbi.1005183) as an extension of the Universal Pacemaker (UPM). The EPM can model non-linear 
epigenetic trait associations directly without transformation of the phenotype of interest[5](https://doi.org/10.1080/15592294.2019.1623634).

## Installation

```shell
pip3 install EpigeneticPacemaker
``` 

## Documentation

[epigeneticpacemaker.readthedocs.io](https://epigeneticpacemaker.readthedocs.io/en/latest/)

## Citations 

1. [Horvath, S. DNA methylation age of human tissues and cell types. Genome Biol. 14, R115 (2013).](https://doi.org/10.1186/gb-2013-14-10-r115)
2. [Orozco, L. D. et al. Epigenome-wide association in adipose tissue from the METSIM cohort. Hum. Mol. Genet. 0, 223495 (2018).](https://doi.org/10.1093/hmg/ddy093)
3. [Snir, S. & Pellegrini, M. An epigenetic pacemaker is detected via a fast conditional expectation maximization algorithm. 10, 695–706 (2018).](https://doi.org/10.1371/journal.pcbi.1005183)
4. [Snir, S., vonHoldt, B. M. & Pellegrini, M. A Statistical Framework to Identify Deviation from Time Linearity in Epigenetic Aging. PLoS Comput. Biol. 12, 1–15 (2016).](https://doi.org/10.2217/epi-2017-0130)
5. [Snir, S., Farrell, C. & Pellegrini, M. Human epigenetic ageing is logarithmic with time across the entire lifespan. Epigenetics (2019). doi:10.1080/15592294.2019.1623634](https://doi.org/10.1080/15592294.2019.1623634)
