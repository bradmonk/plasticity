PWD=$(shell pwd)

help:
	@echo 'Makefile for Cython code                                                    '
	@echo '                                                                            '
	@echo 'Usage:                                                                      '
	@echo '   make _cython_interface.c       builds cython_interface C file            '
	@echo '   make libadvanceonestep.so      builds cython_interface shared object file'
	@echo '   make advance_one_step.mexa64   builds MEX file which uses Cython code    '
	@echo '   make test_mex                  tests the newly built MEX file            '

_cython_interface.c: _cython_interface.pyx
	cython -o _cython_interface.c _cython_interface.pyx

libadvanceonestep.so: _cython_interface.c
	gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
	-I/usr/include/python2.7 -o libadvanceonestep.so _cython_interface.c \
	-lpython2.7

advance_one_step.mexa64: libadvanceonestep.so
	mex -v -I/usr/include/python2.7 -L/home/dhermes/plasticity/PyDiffusion \
	CFLAGS="-fPIC -Wl,-rpath=/home/dhermes/plasticity/PyDiffusion -Wall" \
	-ladvanceonestep -lpython2.7 advance_one_step.c

test_mex:
	LD_LIBRARY_PATH=$(PWD) matlab -nojvm -nodisplay -nosplash -r 'test_mex'

.PHONY: help build_mex test_mex
