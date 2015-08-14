# make C spikefit/gluttonfit library and its MEX interface.
# Notes: sf.c is MEX gateway to spikefit.o libraries.
#        gf.c is MEX gateway to gluttonfit.o libraries.
# todo: make openmp a make task.
# Barnett 2/13/15. glutton 5/4/5, linking spikefit.h 5/19/15
# split off mwrap away from the user's default make 6/18/15
# Tidy 8/14/15

# If you want to use Mwrap, point to it here:
MWRAP = mwrap

.c.o:
	gcc -fPIC -fopenmp -O3 -c $<
#	gcc -fPIC -O3 -c $<

default: greedy glutton

mwrap: sf.c gf.c greedy glutton

sf.c: spikefit.mw
# note -mb here makes separate .m file for each @ function in the .mw ...
	$(MWRAP) -list -mex sf -mb spikefit.mw
	$(MWRAP) -mex sf -c sf.c spikefit.mw

greedy: sf.c spikefit.o spikefit.h
	mex sf.c spikefit.o -lm -lgomp
#	mex sf.c spikefit.o -lm

gf.c: gluttonfit.mw
	$(MWRAP) -list -mex gf -mb gluttonfit.mw
	$(MWRAP) -mex gf -c gf.c gluttonfit.mw

glutton: gf.c gluttonfit.o spikefit.o spikefit.h
	mex gf.c gluttonfit.o spikefit.o -lm -lgomp

clean:
	rm -f *.o *.mex* *~
# note, doesn't clean out the .c and .m mwrap-generated files
