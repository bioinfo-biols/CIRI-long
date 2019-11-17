## CIRI-long

[![Build Status](https://travis-ci.com/Kevinzjy/CIRI-long.svg?token=sq9vq2uzixNezTJhts8Z&branch=master)](https://travis-ci.com/Kevinzjy/CIRI-long)
[![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](https://github.com/Kevinzjy/CIRI-long/blob/master/LICENSE)

Circular RNA Identification for Nanopore Sequencing 

#### Dependency

- **gcc 4.8+ or clang 3.4+ and cmake 3.2+ is needed**
- **Only python3 is supported**
- CIRI-long requires pysam lib, which need executable binary and header of zlib, bzip2, xz, 
please refer to documentation of [pysam](https://pysam.readthedocs.io/en/latest/installation.html) for installation instructions
- all python dependencies are listed in `requirements.txt`


#### Installation

You can simply clone the whole repository, then use `make` to start a complete installation 

```bash
git clone https://github.com/Kevinzjy/CIRI-long.git CIRI-long
cd CIRI-long

# Install submodule
git submodule init 
git submodule update

# Install CIRI-long 
make
```

**Note**: for expert users, the installation under virtualenv is highly recommended

```bash
git clone https://github.com/Kevinzjy/CIRI-long.git CIRI-long
cd CIRI-long

# Install submodule
git submodule init 
git submodule update

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

