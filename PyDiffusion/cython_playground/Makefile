PWD=$(shell pwd)

help:
	@echo 'Makefile for Cython code                                       '
	@echo '                                                               '
	@echo 'Usage:                                                         '
	@echo '   make foo.c        builds cython_interface C file            '
	@echo '   make libfoo.so    builds cython_interface shared object file'
	@echo '   make test         builds test binary                        '
	@echo '   make run_test     runs test binary                          '
	@echo '   make clean        cleans all generated outputs              '

foo.c: foo.pyx
	cython -o foo.c foo.pyx

libfoo.so: foo.c
	gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
	-I/usr/include/python2.7 -o libfoo.so foo.c \
	-lpython2.7

test: libfoo.so main.c
	gcc -I/usr/include/python2.7 \
	-L$(PWD) \
	-Wl,-rpath=$(PWD) -Wall \
	-o test main.c -lfoo -lpython2.7

run_test: test
	PYTHONPATH='$${PYTHONPATH}:$(PWD)' ./test

clean:
	rm -f bar.pyc foo.c foo.h libfoo.so test

.PHONY: help run_test clean
