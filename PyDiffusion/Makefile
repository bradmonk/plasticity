PWD=$(shell pwd)

# On OS X, we need to do a little extra work for the MEX linker.
UNAME=$(shell uname -s)
ifeq ($(UNAME), Darwin)
	MEX_EXTRA="-L/System/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/config"
endif

help:
	@echo 'Makefile for Cython code                                                    '
	@echo '                                                                            '
	@echo 'Primary Targets:                                                            '
	@echo '   make test_mex                  tests the newly built MEX file            '
	@echo '   make matlab                    runs matlab with the MEX/Python paths set '
	@echo '                                                                            '
	@echo 'Other Build Targets:                                                        '
	@echo '   make _cython_interface.c       builds cython_interface C file            '
	@echo '   make libadvanceonestep.so      builds cython_interface shared object file'
	@echo '   make advance_one_step.mexa64   builds MEX file which uses Cython code    '
	@echo '   make clean                     cleans all generated outputs              '

_cython_interface.c: _cython_interface.pyx
	cython -o _cython_interface.c _cython_interface.pyx

libadvanceonestep.so: _cython_interface.c
	gcc -shared -pthread -fPIC -fwrapv -O2 -Wall -fno-strict-aliasing \
	-I/usr/include/python2.7 -o libadvanceonestep.so _cython_interface.c \
	-lpython2.7

advance_one_step.mexa64: advance_one_step.c libadvanceonestep.so
	mex -L$(PWD) $(MEX_EXTRA) -I/usr/include/python2.7 \
	-ladvanceonestep -lpython2.7 -ldl advance_one_step.c

test_mex: advance_one_step.mexa64
	PYTHONPATH="$${PYTHONPATH}:$(PWD)" LD_LIBRARY_PATH=$(PWD) matlab \
	-nojvm -nodisplay -nosplash -r 'test_mex; exit;'

matlab: foo_mex.mexa64
	PYTHONPATH="$${PYTHONPATH}:$(PWD)" LD_LIBRARY_PATH=$(PWD) matlab \
	-nojvm -nodisplay -nosplash

clean:
	rm -f _cython_interface.c _cython_interface.h libadvanceonestep.so \
	 advance_one_step.mexa64

.PHONY: help test_mex matlab clean
