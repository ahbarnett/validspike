% make synth data to make original paper Fig 7b EJ fit.
% Barnett 1/29/16

if 0   % extract wf from the fit for Fig 7b...
  tic; load data_valid/tseries_thresh100K10keep10.mat; toc
  %Elapsed time is 78.881298 seconds.    ceph problem

  wf.dinfo = rmfield(wf.dinfo,'b');
  wf.dinfo = rmfield(wf.dinfo,'sig');
  K = size(wf.W,3);
  wf.freqs = histc(o.L,1:K)  / d.T;   % overwrite with true rates (from fit)
  save data_valid/wf_for_synth_matching_fig7b.mat wf
end

