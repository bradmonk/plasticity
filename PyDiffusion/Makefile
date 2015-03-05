PWD=$(shell pwd)

help:
	@echo 'Makefile for Cython code                                             '
	@echo '                                                                     '
	@echo 'Usage:                                                               '
	@echo '   make cython_interface   builds cython_interface shared object file'
	@echo '   make build_mex          builds MEX file which uses Cython code    '

cython_interface:
	cython -o _cython_interface.c _cython_interface.pyx
	gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
	-I/usr/include/python2.7 -o _cython_interface.so _cython_interface.c \
	-lpython2.7

build_mex:
	cp _cython_interface.so libadvanceonestep.so
	mex -I/usr/include/python2.7 -L/home/dhermes/plasticity/PyDiffusion \
	CFLAGS="-fPIC -Wl,-rpath=/home/dhermes/plasticity/PyDiffusion -Wall" \
	-ladvanceonestep -lpython2.7 advance_one_step.c

.PHONY: help cython_interface build_mex
