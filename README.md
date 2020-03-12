## CIRI-long

[![Build Status](https://travis-ci.com/Kevinzjy/CIRI-long.svg?token=sq9vq2uzixNezTJhts8Z&branch=master)](https://travis-ci.com/Kevinzjy/CIRI-long)
[![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](https://github.com/Kevinzjy/CIRI-long/blob/master/LICENSE)

Circular RNA Identification for Nanopore Sequencing Data

### Author ###

Authors: Jinyang Zhang(zhangjinyang@biols.ac.cn), Fangqing Zhao(zfdq@biols.ac.cn)

Maintainer: Jinyang Zhang

### Release Notes ###

- version 0.4.0: add support for fetching circRNAs without full-length structure 
- version 0.3.0: stringent filter for repeat segments identification
- version 0.2.0: use minimap / bwapy for bsj detection
- Version 0.1.0: basic function of CIRI-long

### License ### 

The code is released under the MIT License. See the `LICENSE` file for more detail

#### Dependency

- **`gcc 4.8+` or `clang 3.4+` and `cmake 3.2+` is needed**
- **Only python3 is supported**
- CIRI-long requires pysam lib, which need executable binary and header of zlib, bzip2, xz, 
please refer to documentation of [pysam](https://pysam.readthedocs.io/en/latest/installation.html) for installation instructions
- all python dependencies are listed in `requirements.txt`
- `samtools` version 1.9 or higher

#### Installation

You can simply clone the whole repository, then use `make` to start a complete installation 

```bash
# Clone CIRI-long and submodules
git clone --recursive https://github.com/Kevinzjy/CIRI-long.git CIRI-long
cd CIRI-long

# Install CIRI-long 
make

# Test CIRI-long
make test
```

**Note**: for expert users, the installation under virtualenv is highly recommended

```bash
git clone --recursive https://github.com/Kevinzjy/CIRI-long.git CIRI-long
cd CIRI-long

# Create virtual environment
python3 -m venv venv

# Activate virtualenv
source ./venv/bin/activate

# Install CIRI-long
make

# Test for installation
make test

# Deactivate
deactivate
```

