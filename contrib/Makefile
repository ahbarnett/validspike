# Makefile for contributed packages
# Barnett 6/4/15. mwrap removed 6/18/15. clean 8/13/15

default:
	mex pdistmex.cpp
	(cd dbscan_patwary; make)
clean:
	rm -f pdistmex.mex*
	(cd dbscan_patwary; make clean)
