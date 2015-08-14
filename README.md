# validspike: spike sorting and validation tools in MATLAB

This is a package for spike sorting of extracellular multichannel neuronal recordings, computing validation metrics, and generating synthetic datasets.

  Alex Barnett and Jeremy Magland, November 2014 - August 2015, SCDA, Simons Foundation
  ahb@math.dartmouth.edu

This accompanies the paper, "Validation of neural spike sorting algorithms without ground-truth information", Alex H. Barnett, Jeremy F. Magland, and Leslie F. Greengard (2015)

## Dependencies

This has been developed on a linux environment running MATLAB R2014a without any toolboxes.
The `mex` compiler is used. Also needed: GNU make, GCC, openmp.
(It has not been tested on Mac OSX; here you may need to adjust compiler options and the location of the `mex` compiler command.)
Optional dependencies: [spikespy](https://github.com/magland/spikespy) for time-series viewing, and
[MWrap](http://www.cs.cornell.edu/~bindel/sw/mwrap/) if you want to generate your own MEX interfaces.

## Setting up

1. If you don't have a MATLAB `startup.m` file, create one as explained [here](http://www.mathworks.com/help/matlab/ref/startup.html). If you want to be able to use the optional `spikespy` viewer in the example scripts, edit the path to it in `vs_startup.m`. Now include the following line in your `startup.m`, replacing `VSHOME` by the full path for this `validspike` package:

        run VSHOME/vs_startup

2. In the top level (`validspike`) directory, run `make` to compile various C/MEX codes and contributed codes. You may need to adjust compiler flags in the file `mexopts.sh` in MATLAB's installation. For instance, on an ubuntu 12.04 system running MATLAB R2012a I needed to add the flag `--openmp` to `LDFLAGS` in the `glnxa64)` case.

3. From MATLAB, run `testall` which should take 1-2 minutes and complete without error. It will generate synthetic data files of 100 MB in size in the `data/` directory from the example waveforms there. These synthetic files are then used in the self-tests of various functions.

See the data load function `examples/loaddata.m` to see ideas for loading your own datasets. Currently this function is set up to point to various datasets in `data_external/` which are not distributed.

## Directories in this package

The directories are as follows, organized by task. We also list some but not all of the main routines in each directory.
In MATLAB you can also use eg `help stageA` to see a list of routines in directory `stageA`.
Most routines can be self-tested by running without arguments. See `testall.m` for a list of all tests.
The acronym EC means extra-cellular, the type of time-series recording that spike sorting
applies to.

####synthesis

Codes for generating artificial datasets given a set of waveform shapes

`synth_Poissonspiketrain` make synthetic time-series with noise  
`synth_singleeventclips`  make synthetic single-spike events (clips) with noise  
`loaddemodata`  loads to memory a default synthetic extracellular times-series  
`loaddemoclips` loads to memory a default synthetic set of upsampled, aligned clips  


####stageA

Algorithms for filtering, detection of events, basic preprocessing (eg measuring noise)

####stageB

Spike sorting algorithms for clustering, extraction of waveforms, clip-based sorting

####stageC_fitlib

Time-series fitting algorithms (greedy, etc) that assume a known waveform set

####validation

Algorithms for validating spike-sorting accuracy and stability

`labels_accuracy`    summarize accuracy of a list of labels against "truth" list  
`times_labels_accuracy`    summarize accuracy of a list of times and labels against "truth" such lists  
`labels_similarity`  report confusion matrix and accuracy stats between two lists of labels  
`getbestshuffling`  find best label permutation between two integer lists  
`times_labels_confusion_matrix`  best-permutated confusion matrix between two lists of times and labels  

####utilities

Assorted low-level functions  

####visualization

Note: many better visualization features are available in the `spikespy' package.  
`plot_spike_shapes`  plots a set of multi-channel waveforms  
`listenraw`  play raw EC data as audio & write out .WAV file  
`viewraw`  plot raw EC data; replaced by `spikespy` package

####examples

These are driver scripts that can be called without arguments, and run complete experiments

`driver_clips`  demonstrate clip-based spike sorting on synthetic clips data  
`driver_timeseries`  demonstrate spike sorting on synthetic time series data  
`loaddata` an example of loading EC timeseries object from various (not supplied) datasets  

####contrib

3rd-party codes and utilites from other sources

## Other directories

`data/`  example waveforms, and where around 100 MB of synthetic data is placed  
`paper_fig_drivers/` figure-generating codes for the above preprint and other research  
`data_external/` a data directory used by us, that is accessed by some advanced examples and figure-generating codes  
