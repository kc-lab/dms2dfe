#dms2dfe

Overview
--------

dms2dfe (**D**eep **M**utational **S**canning (DMS) data to **D**istribution of **F**itness **E**ffects (DFE)) is designed for analyzing Deep Mutational Scaning [1] data in terms of Distribution of Fitness Effects.
In cases where biology do not involve actual fitness, fitness can be regarded as a proxy for preferential enrichment.
Here 'Fitness' and 'selection' here basically are conceptually similar to 'fold changes' and 'up/down regulations' in case of differential expression estimations respectively.

Basic usage
-----------

##Installation

	git clone https://github.com/kc-lab/dms2dfe.git
	cd dms2dfe
	pip install .

##Usage

    import dms2dfe
    dms2dfe.pipeline.main(`path/to/project_directory`)

##Documentation

Documentation of dms2dfe is available at kc-lab.github.io/dms2dfe .

Citations
---------

[1] Fowler, D.M., and S. Fields. 2014. Deep mutational scanning: a new style of protein science. Nature methods. 11: 801–7.
