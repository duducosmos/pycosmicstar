#!/usr/bin/env python
# -*- coding: utf-8 -*-


# metadata
""" The Cosmic Start Formation Rate is a software to study the star
formation history for diferents cosmological models """
__version__ = " 0.1 "
__license__ = " GNU GENERAL PUBLIC LICENSE V3 "
__author__ = " Eduardo dos Santo Pereira "
__email__ = " pereira.somoza@mail.com "
__url__ = " www.cosmicstarformation.com "
__date__ = " 2013-12-02T15:58:52 "
__prj__ = " thecosmicstarformationrate "


# imports
import os,glob
from setuptools import setup, find_packages
from numpy.distutils.core import setup, Extension

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


Extensions = [Extension(name="lcdmlib",
                        sources=["cosmicstarformation/lcdmlib.f"],)]


setup(
    name="cosmicstarformation",
    version="0.1",
    url="www.cosmicstarformation.com",
    download_url="www.cosmicstarformation.com",
    license="GNU GENERAL PUBLIC LICENSE V3",
    author="Eduardo dos Santo Pereira",
    author_email="pereira.somoza@mail.com",
    maintainer="Eduardo dos Santo Pereira",
    maintainer_email="pereira.somoza@mail.com",
    description="""The Cosmic Start Formation Rate is a software to study the
    star formation history for diferents cosmological models""",
    packages=find_packages(),
    package_data={'': ['*.so']},
    #data_files=datafiles,
    install_requires=['numpy', 'scipy'],
    long_description=read('README'),
    ext_modules=Extensions,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: GNU V3",
    ],


)