# CIRI-long
Circular RNA Identification for Nanopore Sequencing 

```bash
git clone --recursive https://github.com/rvaser/spoa spoa
cd spoa
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cmake -DCMAKE_INSTALL_PREFIX:PATH=/usr . && make all install

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
```