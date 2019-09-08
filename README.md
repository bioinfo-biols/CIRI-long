# CIRI-long
Circular RNA Identification for Nanopore Sequencing 

#### Installation

You can simply clone the whole repository, then use `make` to start a complete installation 

```bash
git clone --recursive https://github.com/Kevinzjy/CIRI-long.git CIRI-long
cd CIRI-long
make
```

**Note**: for expert users, the installation under virtualenv is highly recommended

```bash
git clone --recursive https://github.com/Kevinzjy/CIRI-long.git CIRI-long
cd CIRI-long

# Create virtual environment
python -m venv venv

# Activate virtualenv
source ./venv/bin/activate

# Install CIRI-long
make

# Deactivate
deactivate
```