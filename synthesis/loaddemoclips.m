function [X fac d] = loaddemoclips
% LOADDEMOCLIPS - load (and if needed generate) default upsampled aligned clips
%
% [X fac d] = loaddemoclips loads the default upsampled aligned clips data
%   from data/ directory (it also auto-generates it if absent).
%
% Auto-generation of such clips is done by downsampling from the default
%   waveforms, adding noise, then by upsampling and aligning.
%
% Output:
%  X - M*T*Nc double array of clip data
%  fac - upsampling factor used
%  d - info about the supposed originating EC dataset (d.samplefreq, etc)

% Barnett 6/10/15

if nargout<1, test_loaddemoclips; return; end

h = fileparts(mfilename('fullpath')); % current dir
synfile = fullfile(h,'../data/clips_default_synth.mat');
if ~exist(synfile,'file')                % auto-generate it
  fprintf('generating %s ...\n',synfile)
  wf = loaddefaultwaveforms;
  d = wf.d;
  [d.M,~,K] = size(wf.W);
  Nt = 40;           % 2 ms per synthetic sample-rate clip
  pops = 1e3 * ones(1,K);   % how many clips of each type
  eta = 40;    % was 25 roughly for this dataset
  noi = setup_noisemodel(d,Nt,eta,0.0002);  % choose noise: std dev eta, tau_corr
  o.ampl = 0.1;
  tpad = (Nt-5)/2;   % add small jitter over 5 samples
  [X l] = synth_singlespikeclips(wf, Nt, pops, noi, tpad, 0, o);
  % X is now fake detected clips at the sample rate and with jitter. Now u & a:
  fac = 3;   % choose upsampling (needn't be same as for the originating wf)
  X = upsample(X, fac);
  X = alignspikes(X, fac);
  d.name = ['synth ua clips from ' d.name];
  save(synfile,'X','fac','d');
  gndfile = fullfile(h,'../data/clips_default_synth_groundtruth.mat');
  save(gndfile,'l');   % labels only
else
  load(synfile);
end
%%%%
  
function test_loaddemoclips
tic; [X fac d] = loaddemoclips; fprintf('load time = %.3g s\n',toc)
'demo clips have X size:', size(X),fac,d
plot_spike_shapes(X(:,:,1:100),sprintf('1st 100 clips from %s',d.name));
