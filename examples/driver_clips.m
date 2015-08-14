% example script for spike-sorting a set of upsampled & aligned clips
% Barnett 6/10/15

clear
[X fac d] = loaddemoclips;    % get some upsampled+aligned clips in X
o.cmethod = 'k++'; o.K = 7;   % choose clustering method (this one needs K)
[l W] = spikesort_clips(X,o); % assumes upsampled+aligned

overview_sorted_clips(X,l);   % show some clips grouped by sorted label
plot_spike_shapes(W,'waveforms found');

% assess accuracy, since we have luxury of ground-truth...
g = load('clips_default_synth_groundtruth.mat');
labels_accuracy(g.l,l);  % accuracy is high since it was easy synthetic data
