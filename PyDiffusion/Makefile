help:
	@echo 'Makefile for Cython code                                             '
	@echo '                                                                     '
	@echo 'Usage:                                                               '
	@echo '   make cython_interface   builds cython_interface shared object file'

cython_interface:
	cython -o _cython_interface.c _cython_interface.pyx
	gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
	-I/usr/include/python2.7 -o _cython_interface.so _cython_interface.c

.PHONY: cython_interface help
