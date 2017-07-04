
include make.inc

all: libdfftpack testdfftpack

libdfftpack:
	mkdir -p ./lib
	mkdir -p ./objs
	( cd ./src; $(MAKE) clean; $(MAKE) )

testdfftpack:
	( cd ./test; $(MAKE) clean; $(MAKE) )

clean:
	( cd ./src; $(MAKE) clean; cd ../test; $(MAKE) clean )
