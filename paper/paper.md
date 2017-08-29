---
title: 'dms2dfe: Comprehensive Workflow for Analysis of Deep Mutational Scanning Data'
authors:
- affiliation: 1
  name: Rohan Dandage
  orcid: 0000-0002-6421-2067
- affiliation: 1
  name: Kausik Chakraborty
  orcid: 0000-0001-6000-8379
affiliations:
- index: 1
  name: CSIR-Institute of Genomics and Integrative Biology, New Delhi, India.
date: "04 July 2017"
output: pdf_document
bibliography: paper.bib
tags:
- deep sequencing
- mutations
---

# Summary

`dms2dfe` is a python package that integrates steps in the analysis of Deep Mutational Scanning [@Fowler2014a] data. Using this end-to-end workflow, users can implement various processing methods and downstream applications in the deep sequencing based bulk mutational scanning experiments.

Recently, owing to evolution of sequencing and phenotyping technologies, large scale genotype to phenotype data is increasingly being generated. Along this line of research, Deep Mutational Scanning method allows comprehensive assessment of all the substitutions of a given gene. In the analysis of Deep Mutational Scanning data, `dms2dfe` allows end-to-end workflow that addresses issue of noise control and offers variety of crucial downstream analyses. In downstream analyses, `dms2dfe` workflow provides estimation of distributions of effect sizes, identification of potential molecular constraints and generation of data-rich visualizations [@Dandage2016]. While a number of tools have been developed for analysis of DMS data [@Fowler2011;@Bloom2015;@Rubin2016b], `dms2dfe` workflow integrates various state-of-art analysis tools such as Trimmomatic [@Bolger2014], Bowtie [@Langmead2012], samtools [@Li2009] and DESeq2 [@Love2014] for an end-to-end pairwise analysis of samples. Users familiar with such widely used tools can opt for `dms2dfe` workflow.

As an input for the workflow, deep sequencing data (whether unaligned or aligned) or list of genotypic varants can be provided. For data structure, `dms2dfe` uses DataFrames from robust Pandas library [@McKinney2010]. For a demonstration purpose, sample datasets from various studies [@Firnberg2014;@Olson2014;@Melnikov2014a] are available here ^[<https://github.com/rraadd88/ms_datasets>].

Source code and issue tracker is available in `dms2dfe`'s GitHub repository ^[<https://github.com/kc-lab/dms2dfe>]. Documentarion and API ^[<https://kc-lab.github.io/dms2dfe>] are generated using Sphinx ^[<http://www.sphinx-doc.org>].

# References