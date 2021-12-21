.PHONY: help all lib prepare install test clean

all: lib prepare install

help:
	@echo "make all"
	@echo "       install CIRI-long and dependencies into current python environment"
	@echo "make lib"
	@echo "       prepare development environment, use only once"
	@echo "make prepare"
	@echo "       install python requirements"
	@echo "make install"
	@echo "       install CIRI-long"
	@echo "make test"
	@echo "       unit test, run after installation"
	@echo "make clean"
	@echo "       clean python cache files"

lib:
	cd libs/striped_smith_waterman && make libssw.so

prepare:
	cat requirements.txt | xargs -n 1 pip3 install

install:
	python3 setup.py install

test:
	python3 setup.py build_ext --inplace
	python3 setup.py test

clean:
	python3 setup.py clean --all
