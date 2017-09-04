---
title: 'dms2dfe: Comprehensive Workflow for Analysis of Deep Mutational Scanning Data'
tags:
  - deep sequencing
  - mutations
authors:
  - name: Rohan Dandage
    orcid: 0000-0002-6421-2067
    affiliation: 1
  - name: Kausik Chakraborty
    orcid: 0000-0001-6000-8379
    affiliation: 1
affiliations:
  - index: 1
    name: CSIR-Institute of Genomics and Integrative Biology, New Delhi, India.
date: "04 July 2017"
output: pdf_document
bibliography: paper.bib
---
<!--
pandoc --bibliography=paper.bib --template=latex.template paper.md -o paper.pdf
-->

# Summary

`dms2dfe` is an integrative analysis workflow designed for end-to-end analysis of Deep Mutational Scaning [@Fowler2014a] data. Using this workflow, users can implement various processing methods and downstream applications for pairwise enrichment analysis of ultra-deep sequencing data.

Recently, owing to the evolution of sequencing and phenotyping technologies, large scale genotype to phenotype data is increasingly being generated. Along this line of research, Deep Mutational Scanning method allows comprehensive assessment of all the possible amino acid substitutions of an entire gene or part of a gene. In the analysis of Deep Mutational Scanning data, `dms2dfe` addresses crucial issue of noise control using widely used DESeq2 [@Love2014] workflow and offers variety of downstream analyses to contextualize the results. In downstream analyses, `dms2dfe` workflow provides identification of potential molecular constraints, comparative analysis across different experimental conditions and generation of data-rich visualizations [@Dandage2016]. While a number of tools have been developed for analysis of DMS data [@Fowler2011;@Bloom2015;@Rubin2016b], users familiar with commonly used state-of-art genomics tools such as Trimmomatic [@Bolger2014], Bowtie [@Langmead2012], samtools [@Li2009] and DESeq2 [@Love2014] can opt for `dms2dfe` workflow for analysis of preferential enrichments. Note that `dms2dfe` workflow is designed exclusively for experimental designs in which there is a need of pairwise analysis of samples eg. before and after selection.

As an input for the workflow, deep sequencing data (whether unaligned or aligned) or list of genotypic varants can be provided.  For a demonstration purpose, sample datasets from various studies [@Firnberg2014;@Olson2014;@Melnikov2014a] are available here ^[<https://github.com/rraadd88/ms_datasets>]. `dms2dfe` uses DataFrames from robust Pandas library [@McKinney2010] for processing all the tabular data.

Source code and issue tracker is available in `dms2dfe`'s GitHub repository ^[<https://github.com/kc-lab/dms2dfe>]. Documentarion and API ^[<https://kc-lab.github.io/dms2dfe>] are generated using Sphinx ^[<http://www.sphinx-doc.org>].

# References