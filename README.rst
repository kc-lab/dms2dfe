``dms2dfe``
===========

|build status| |PyPI version|

Overview
--------

dms2dfe (**D**\ eep **M**\ utational **S**\ canning (DMS) data to
**D**\ istribution of **F**\ itness **E**\ ffects (DFE)) is an
integrative analysis workflow designed for end-to-end enrichment
analysis of Deep Mutational Scaning [1] data. Among many experimental
designs Deep Mutational Scaning asssay, this workflow, currently, is
only applicable for **pairwise comparisons** of samples.

Full documentation
------------------

Latest version: http://kc-lab.github.io/dms2dfe

dms2dfe v1.0.6: http://kc-lab.github.io/dms2dfe/v1.0.6/html/

Quick installation
------------------

To install the package written in python 2.7, simply execute following
command:

::

    pip install dms2dfe

Quick usage
-----------

From bash command line, create a project directory

::

    dms2dfe project_directory

Insert input parameters in the configuration files (.csv) located in
``project_directory/cfg``

Run the analysis,

::

    dms2dfe project_directory

Publication
-----------

| **dms2dfe: Comprehensive Workflow for Analysis of Deep Mutational
  Scanning Data**
| Rohan Dandage, Kausik Chakraborty
| doi: http://dx.doi.org/10.1101/072645

Questions
---------

Please mention them here: https://github.com/kc-lab/dms2dfe/issues .

Footnotes
---------

[1] Fowler, D.M., and S. Fields. 2014. Deep mutational scanning: a new
style of protein science. Nature methods. 11: 801–7.

.. |build status| image:: http://img.shields.io/travis/rraadd88/dms2dfe/master.svg?style=flat
   :target: https://travis-ci.org/rraadd88/dms2dfe
.. |PyPI version| image:: https://badge.fury.io/py/dms2dfe.svg
   :target: https://pypi.python.org/pypi/dms2dfe
