function d=loaddemodata
% LOADDEMODATA - load (and maybe generate) the default EC timeseries dataset
%
% d = loaddemodata loads into the extracellular struct d the default dataset
%  from the data/ directory (it also auto-generates it if absent).
%
% Output struct d contains at least these fields:
%  A             : M (# electrodes) by N (# time samples) extracellular signal
%                  data (arb units)
%  samplefreq    : sampling frequency (Hz)
%  dt            : 1/samplefreq (s)
%  T             : total time (s)
%  electrodelocs : 2*M list of xy-coords of electrodes in order of rows of A
%                  (optional)
%  name          : string description
%
% This defines the EC data struct format
%
% Barnett 6/4/15

if nargout<1, test_loaddemodata; return; end

h = fileparts(mfilename('fullpath')); % current dir
synfile = fullfile(h,'../data/EC_default_synth.mat');
if ~exist(synfile,'file')                % auto-generate it
  fprintf('generating %s ...\n',synfile)
  wf = loaddefaultwaveforms;
  d = wf.d;
  d.dt = 1/d.samplefreq;
  N = round(60*d.samplefreq);            % 1 minute of time
  d.T = N*d.dt;                          % overwrite total time
  noi = setup_noisemodel(d,N,25);        % choose noise std dev eta
  [d.M,~,K] = size(wf.W);
  rates = 20*ones(1,K);                  % spikes/sec for each type; 20 is low
  [d.A p] = synth_Poissonspiketrain(wf,N,rates,noi,2,0); % includes noise
  d.name = 'default synthetic EC data';
  z = [0, exp(2i*pi*(1:6)/6)]; d.electrodelocs = [real(z);imag(z)]; % EJ 7-elec
  save(synfile,'d');
  gndfile = fullfile(h,'../data/EC_default_synth_groundtruth.mat');
  save(gndfile,'p');
else
  load(synfile,'d');
end

%%%%

function test_loaddemodata  % load default data. also tests viewraw
tic; d = loaddemodata; fprintf('load time = %.3g s\n',toc)
viewraw(d);
