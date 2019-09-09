#! /usr/bin/env python
# -*- encoding:utf-8 -*-
import os
import sys

sys.path.append('CIRI')

import codecs
import pathlib

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize
from CIRI.version import __version__


def read(infile):
    return codecs.open(os.path.join(os.path.dirname(__file__), infile)).read()


extensions = [
    Extension(
        "poa",
        sources=["CIRI/poa.pyx"],
        depends=["CIRI/cpoa.h", "CIRI/cpoa.pxd"],
        include_dirs=['./CIRI', 'spoa/include'],
        libraries=['spoa'],
        language="c++",
        # extra_objects=['spoa/build/lib/libspoa.a'],
        library_dirs=['spoa/build/lib'],
        extra_compile_args=['-std=c++11'],
    ),
]

setup(
    name='CIRI-long',
    version=__version__,
    url='https://github.com/Kevinzjy/CIRI-long',
    description='circular RNA identification from Nanopore',
    long_description=read('README.md'),
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
    ext_modules=cythonize(extensions),
    include_package_data=True,
    zip_safe=False,
    install_requires=[
        'argparse==1.2.1', 'Cython==0.29.13', 'mappy==2.17', 'numpy==1.17.0', 'pandas==0.25.0',
        'pysam==0.15.3', 'python-Levenshtein==0.12.0', 'scikit-bio==0.5.5',
        'scikit-learn==0.21.3', 'scipy==1.3.1',
    ],
    test_suite="nose.collector",
    tests_require=['nose'],
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