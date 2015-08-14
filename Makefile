# Makefile for spike-sorting toolbox
# Barnett 6/4/15. clean 8/13/15

default:
	(cd stageC_fitlib; make)
	(cd contrib; make)

clean:
	rm -f data/*default_synth*
	find . -type f -name '*~' -delete
	(cd stageC_fitlib; make clean)
	(cd contrib; make clean)
