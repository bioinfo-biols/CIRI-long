#! /usr/bin/env python
# -*- encoding:utf-8 -*-
import os
import sys

sys.path.append('CIRI')

import codecs


from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize
from CIRI.version import __version__


def read(infile):
    return codecs.open(os.path.join(os.path.dirname(__file__), infile)).read()


setup(
    name='CIRI-long',
    version=__version__,
    url='https://github.com/bioinfo-biols/CIRI-long',
    description='circular RNA identification from Nanopore',
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    author='Jinyang Zhang',
    author_email='zhangjinyang@biols.ac.cn',
    maintainer='Jinyang Zhang',
    maintainer_email='zhangjinyang@biols.ac.cn',
    license='MIT',
    keywords='circRNA',
    packages=find_packages(exclude=['doc', 'tests']),
    entry_points={
        'console_scripts': [
            'CIRI-long=CIRI.main:main',
        ]
    },
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'argparse>=1.2.1', 'Cython>=0.29.13', 'mappy>=2.17', 'numpy>=1.17.0', 'pandas>=0.25.0',
        'pysam>=0.15.3', 'python-Levenshtein>=0.12.0', 'scikit-learn>=0.21.3', 'scipy>=1.3.1',
    ],
    test_suite="nose.collector",
    tests_require=['nose==1.3.7'],
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Intended Audience :: Developers',
        'Operating System :: OS Independent',
    ],
    cmdclass={
        'build_ext': build_ext,
    }
)