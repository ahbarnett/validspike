% run all tests and self-tests for validspike (without checking output)
% Barnett 6/8/15. Changed project name, new conf mat tests 8/14/15

clear; 'please wait around 100 seconds... (apologies for the focus grabs)'
t_testall = tic; save('/tmp/testall.mat'); % save the start time info from clears

%%%%%%%%%%%%%%%%%% the tests...

% synthesis
system('rm -f data/*default_synth*');  % make clean the generated default data
loaddemodata
loaddemoclips
setup_noisemodel
noisesample
spikemodel  % older, replaced by spikemod MEX
synth_Poissonspiketrain
synth_singlespikeclips

% stage A
freqfilter
detectevents
signal_moments
% todo: self-tests for empirical*.m

% stage B
alignspikes
upsample
features
cluster
% NB: spikesort_clips is tested by driver_clips

% stage C
test_spikemod
test_fitonesp
test_multifitgreedy
test_synthgreedy
test_fillscore
fit_timeseries
negloglik       % older
fitonespike     % older

% visualization
viewraw
listenraw
plot_spike_shapes
%points_click_waveform  % test waits for user input. obsolete
powerspec
show_crosscorr

% validation
confusion_matrix
bestcolpermconfmat
labels_accuracy
times_labels_confusion_matrix
times_labels_accuracy
spikesetmatch   % legacy, used only to test greedy clips in stageC_fitlib

% utils
% todo: self-tests for utils

% examples
driver_clips
driver_timeseries
driver_clips_stability
driver_timeseries_stability

%%%%%%%%%%%%%%%%%%% clean up...
load('/tmp/testall.mat');  % just to get t_testall
system('rm /tmp/testall.mat');
fprintf('\n...done all validspike tests without crashing in %.3g s; closing figs\n',toc(t_testall))
close all
if exist('spikespy'), spikespy('closeall'); end
