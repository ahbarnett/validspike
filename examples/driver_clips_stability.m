% example script for measuring stability of spike-sorting on a set of upsampled
% & aligned clips.
% Built upon driver_clips.m
% Barnett 6/16/15, tweaked 8/14/15

clear
[X fac d] = loaddemoclips;    % get some ua clips in X
o.K = 7;                      % tell sorter how many spike types to cluster into
o.cmethod='k++'; o.num_trials = 3;  % sorter alg opts
oo.meth='rerun'; oo.num_runs = 10; oo.verb = 2;  % validation opts
[fbar,fsam,info] = eval_stability_clipbased(@spikesort_clips, X, o, oo);
