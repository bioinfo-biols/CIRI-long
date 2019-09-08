.PHONY: all lib CIRI_long test

all: lib CIRI-long

CIRI-long:
	python setup.py install

lib:
	mkdir -p spoa/build
	cd spoa/build && cmake -DCMAKE_BUILD_TYPE=Release .. && make

test:
	python setup.py build_ext --inplace
	mv poa.cpython-36m-x86_64-linux-gnu.so ./CIRI
