PWD=$(shell pwd)

# On OS X, we need to do a little extra work for the MEX linker.
UNAME=$(shell uname -s)
ifeq ($(UNAME), Darwin)
	MEX_EXTRA="-L/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/config"
endif

help:
	@echo 'Makefile for Cython code                                           '
	@echo '                                                                   '
	@echo 'Primary Targets:                                                   '
	@echo '   make run_test         runs test C binary that uses Python       '
	@echo '   make test_mex         tests MEX file that uses Python           '
	@echo '   make matlab           runs matlab with the MEX/Python paths set '
	@echo '                                                                   '
	@echo 'Other Build Targets:                                               '
	@echo '   make foo.c            builds cython_interface C file            '
	@echo '   make libfoo.so        builds cython_interface shared object file'
	@echo '   make test             builds test binary                        '
	@echo '   make foo_mex.mexa64   builds MATLAB mex file                    '
	@echo '   make clean            cleans all generated outputs              '

foo.c: foo.pyx
	cython -o foo.c foo.pyx

libfoo.so: foo.c
	gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
	-I/usr/include/python2.7 -o libfoo.so foo.c -lpython2.7 \

test: libfoo.so main.c
	gcc -I/usr/include/python2.7 \
	-L$(PWD) -Wall \
	-o test main.c -lfoo -lpython2.7 -ldl

run_test: test
	PYTHONPATH="$${PYTHONPATH}:$(PWD)" LD_LIBRARY_PATH=$(PWD) ./test

foo_mex.mexa64: foo_mex.c libfoo.so
	mex -L$(PWD) $(MEX_EXTRA) -I/usr/include/python2.7 -lpython2.7 -lfoo -ldl foo_mex.c

test_mex: foo_mex.mexa64
	PYTHONPATH="$${PYTHONPATH}:$(PWD)" LD_LIBRARY_PATH=$(PWD) matlab -nojvm -nodisplay -nosplash -r 'test_mex; exit;'

matlab: foo_mex.mexa64
	PYTHONPATH="$${PYTHONPATH}:$(PWD)" LD_LIBRARY_PATH=$(PWD) matlab

clean:
	rm -f bar.pyc foo.c foo.h libfoo.so test foo_mex.mexa64

.PHONY: help run_test clean test_mex matlab
