
include ./make.inc

all: libdag

libdag:
	mkdir -p ./lib
	mkdir -p ./objs
	( cd ./src; $(MAKE) )

clean:
	( cd ./src; $(MAKE) clean )
