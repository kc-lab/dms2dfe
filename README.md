#dms2dfe

##Overview

dms2dfe (**D**eep **M**utational **S**canning (DMS) data to **D**istribution of **F**itness **E**ffects (DFE)) is designed for analyzing Deep Mutational Scaning [1] data in terms of Distribution of Fitness Effects.

##Installation

	git clone https://github.com/kc-lab/dms2dfe.git
	cd dms2dfe
	pip install .

##Usage

From python environment,

    import dms2dfe
    dms2dfe.pipeline.main("path/to/project_directory")

From bash command line,

    python path/to/dms2dfe/pipeline.py path/to/project_directory

##Detailed documentation

Available at https://kc-lab.github.io/dms2dfe .

##Citation

**dms2dfe: Comprehensive Workflow for Analysis of Deep Mutational Scanning Data**  
Rohan Dandage, Kausik Chakraborty  
doi: http://dx.doi.org/10.1101/072645  

##Footnotes

[1] Fowler, D.M., and S. Fields. 2014. Deep mutational scanning: a new style of protein science. Nature methods. 11: 801–7.
