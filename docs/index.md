<h1> Epigenetic Pacemaker: A fast conditional expectation maximization algorithm for modeling epigenetic states </h1>

## About
DNA methylation is widely used to model physiological phenotypes, such as 
aging[<sup>1</sup>](https://doi.org/10.1186/gb-2013-14-10-r115) and type II diabetes[<sup>2</sup>](https://doi.org/10.1093/hmg/ddy093). 
The Epigenetic Pacemaker, **EPM**, is an implementation of a fast conditional expectation maximization algorithm that models 
epigenetic states under and evolutionary framework[<sup>3</sup>](https://doi.org/10.2217/epi-2017-0130). The EPM was first introduced by Snir et al.
[<sup>4</sup>](https://doi.org/10.1371/journal.pcbi.1005183) as an extension of the Universal Pacemaker (UPM) model of genome evolution. In contrast to regression 
bases approaches, the EPM does not assume a linear relationship between the epigenetic state and a trait of interest.
As a result the EPM can model non-linear epigenetic trait associations directly without transformation of the phenotype of
interest[<sup>5</sup>](https://doi.org/10.1080/15592294.2019.1623634).
The software implementation of the Epigenetic Pacemaker is described our publication [*The Epigenetic Pacemaker - modeling epigenetic states under an evolutionary framework*](https://academic.oup.com/bioinformatics/article-abstract/doi/10.1093/bioinformatics/btaa585/5861533?redirectedFrom=fulltext)<sup>6</sup>.

## Installation

```shell
pip3 install EpigeneticPacemaker
``` 

<h2> Citations </h2> 

1. [Horvath, S. DNA methylation age of human tissues and cell types. Genome Biol. 14, R115 (2013).](https://doi.org/10.1186/gb-2013-14-10-r115)
2. [Orozco, L. D. et al. Epigenome-wide association in adipose tissue from the METSIM cohort. Hum. Mol. Genet. 0, 223495 (2018).](https://doi.org/10.1093/hmg/ddy093)
3. [Snir, S. & Pellegrini, M. An epigenetic pacemaker is detected via a fast conditional expectation maximization algorithm. 10, 695–706 (2018).](https://doi.org/10.1371/journal.pcbi.1005183)
4. [Snir, S., vonHoldt, B. M. & Pellegrini, M. A Statistical Framework to Identify Deviation from Time Linearity in Epigenetic Aging. PLoS Comput. Biol. 12, 1–15 (2016).](https://doi.org/10.2217/epi-2017-0130)
5. [Snir, S., Farrell, C. & Pellegrini, M. Human epigenetic ageing is logarithmic with time across the entire lifespan. Epigenetics (2019). doi:10.1080/15592294.2019.1623634](https://doi.org/10.1080/15592294.2019.1623634)
6. [Colin Farrell, Sagi Snir, Matteo Pellegrini, The Epigenetic Pacemaker - modeling epigenetic states under an evolutionary framework, Bioinformatics, , btaa585, https://doi.org/10.1093/bioinformatics/btaa585](https://academic.oup.com/bioinformatics/article-abstract/doi/10.1093/bioinformatics/btaa585/5861533?redirectedFrom=fulltext)



