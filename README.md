# `dms2dfe`

[![build status](
  http://img.shields.io/travis/rraadd88/dms2dfe/master.svg?style=flat)](
 https://travis-ci.org/rraadd88/dms2dfe) [![PyPI version](https://badge.fury.io/py/dms2dfe.svg)](https://pypi.python.org/pypi/dms2dfe)

## Overview

dms2dfe (**D**eep **M**utational **S**canning (DMS) data to **D**istribution of **F**itness **E**ffects (DFE)) is designed for analyzing Deep Mutational Scaning [1] data in terms of Distribution of Fitness Effects.

## Full documentation

Version 1.x.x: http://kc-lab.github.io/dms2dfe/v1.0.0/html/  

## Quick installation

Requires python 2.7 environment. Tested on Ubuntu 12.04. To install download the latest release from https://github.com/kc-lab/dms2dfe/releases and execute following commands.  

    cd dms2dfe
    pip install .

## Quick usage

From bash command line, create a project directory

    dms2dfe path/to/project_directory

Insert input parameters in the configuration files (.csv) located in `project_directory/cfg`   

Run the analysis,

    dms2dfe path/to/project_directory


## Publication

**dms2dfe: Comprehensive Workflow for Analysis of Deep Mutational Scanning Data**  
Rohan Dandage, Kausik Chakraborty  
doi: http://dx.doi.org/10.1101/072645  

## Questions

Please mention them here: https://github.com/kc-lab/dms2dfe/issues .

## Footnotes

[1] Fowler, D.M., and S. Fields. 2014. Deep mutational scanning: a new style of protein science. Nature methods. 11: 801–7.
