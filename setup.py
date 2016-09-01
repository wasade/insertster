#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2016--, insertster development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

version = '0.0.1-dev'

classes = """
    Development Status :: 4 - Beta
    License :: OSI Approved :: BSD License
    Topic :: Software Development :: Libraries
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3 :: Only
    Programming Language :: Python :: 3.4
    Programming Language :: Python :: 3.5
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""
classifiers = [s.strip() for s in classes.split('\n') if s]

description = ('Insert sequences into a phylogenetic tree')

with open('README.rst') as f:
    long_description = f.read()

setup(name='insertster',
      version=version,
      license='BSD',
      description=description,
      long_description=long_description,
      author="insertster development team",
      author_email="wasade@gmail.com",
      maintainer="insertster development team",
      maintainer_email="wasade@gmail.com",
      packages=find_packages(),
      install_requires=[
          'scikit-bio >= 0.5.0, < 0.6.0'
      ],
      classifiers=classifiers,
      package_data={
          'insertster.tests': ['data/*'],
          }
      )
