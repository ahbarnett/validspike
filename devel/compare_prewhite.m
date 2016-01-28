% prewhitening tests, on Buzsaki, for stageA/channelprewhiten.m
% Adapted from demo_spikesort_timeseries_buzsaki.m
% Barnett 1/27/16.

clear;
dd = loaddata('b'); fprintf('loaded\n')   % 2.5 min, 10 channel
dd.A = freqfilter(dd.A,dd.samplefreq,300,[]);
%spikespy({d.A,'before white'});

for run=1:3
  wo = []; wo.verb = 1; wo.rownorm = 1;  % always a good idea
  if run==1, wo.meth = 'u'; thresh = Inf;  % jeremy
  elseif run==2, wo.meth = 's'; thresh = Inf;  % alex variant mixing chans
  elseif run==3, wo.meth = 'c'; thresh = [];     % alex from 2015
  end
  d = channelprewhiten(dd,thresh,wo);   % what we're testing
  % normalize each channel (weirdly needed for meth='s' even if Q row-normed):
  for m=1:size(d.A,1), d.A(m,:) = d.A(m,:)/norm(d.A(m,:))*norm(dd.A(1,:)); end
  spikespy({d.A,sprintf('after white: wo.meth=%s',wo.meth)});
  %max(d.A(:))
  o = []; o.verb = 0;         % try clustering (no fitting)...
  o.cmethod='k++'; o.K = 7; o.Kkeep = o.K; o.thresh=120;
  [wf L z] = waveforms_from_timeseries(d.A,d.samplefreq,o);
  plot_labeled_pts(z,L); title(sprintf('z: wo.meth=%s',wo.meth))
  plot_spike_shapes(wf.W,sprintf('W final: wo.meth=%s',wo.meth));
end

% results: 'u' and 'c' v similar, so switch to u since preserves channel idents.
% using noise clips to est S, vs whole data, is very similar.
% no prewhitening is crap.
% 's' is worse, since need to normalize channels afterwards.
