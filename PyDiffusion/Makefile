help:
	@echo 'Makefile for Cython code                               '
	@echo '                                                       '
	@echo 'Usage:                                                 '
	@echo '   make speedup       builds speedup shared object file'

speedup:
	cython -o _speedup.c _speedup.pyx
	gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing -I/usr/include/python2.7 -o _speedup.so _speedup.c

.PHONY: speedup help
