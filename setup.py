#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rohan@igib.in,rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.    


"""
========
setup.py
========

installs dms2dfe

USAGE :
python setup.py install

Or for local installation:

python setup.py install --prefix=/your/local/dir

"""

import sys
try:
    from setuptools import setup, find_packages, Extension 
except ImportError:
    from distutils.core import setup, find_packages, Extension


with open('requirements.txt') as f:
    required = f.read().splitlines()

if (sys.version_info[0], sys.version_info[1]) != (2, 7):
    raise RuntimeError('Python 2.7 required ')
               
# main setup command
setup(
name='dms2dfe',
author='Rohan Dandage',
author_email='rraadd_8@hotmail.com,rohan@igib.in',
version='1.0.6.2',
url='https://github.com/rraadd88/dms2dfe',
download_url='https://github.com/rraadd88/dms2dfe/archive/master.zip',
description='Pipeline to analyse Deep Mutational Scanning (DMS) experiments in terms of Distribution of Fitness Effects (DFE)',
long_description='https://github.com/rraadd88/dms2dfe/README.md',
license='General Public License v. 3',
install_requires=required,
platforms='Linux (Tested on Debian 7 (Ubuntu 12.04))',
keywords=['bioinformatics','Deep sequencing','molecular evolution'],
packages=find_packages(),
package_data={'': ['dms2dfe/tmp','dms2dfe/cfg']},
include_package_data=True,
entry_points={
    'console_scripts': ['dms2dfe = dms2dfe.pipeline:main',],
    },
)


## JUNK ##
##make tags
#git tag -a v$(python setup.py --version) -m "http://kc-lab.github.io/dms2dfe/v1.0.0/html/"
#git push --tags
##rewrite old tags
#git tag new old
#git tag -d old;git push origin :old
#git push --tags
#undo last commit
###git reset --hard HEAD~
## ignore tmp files
# git update-index --assume-unchanged FILE_NAME
# cd docs
# make html
# cd ../
# rsync -avv docs/_build/html docs/latest
# rsync -avv docs/_build/html docs/v1.0.0

# cd ../kc_lab_dms2dfe
# git fetch upstream
# git rebase upstream/master
# git push
# git fetch upstream --tags
# git push --tags
# cd ../dms2dfe

# cd ../kc_lab_io
# rsync -avv ../dms2dfe/docs/latest dms2dfe
# rsync -avv ../dms2dfe/docs/stable dms2dfe
# rsync -avv ../dms2dfe/docs/v1.0.0 dms2dfe
# git add --all
# git commit -m "update"
# git push
# cd ../dms2dfe

## JUNK ##