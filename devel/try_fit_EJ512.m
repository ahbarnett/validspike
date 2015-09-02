% try fitting JE 512 using JFM's waveforms
% Barnett 8/28/15

clear
wf.W = readmda('~/magland/firetrack/testdata/waveforms_first_5e5_points.mda');
%[M T K] = size(wf.W);
wf.fac = 1;   % since no upsampling by JFM

tic; raw = load('~/ss_datasets/EJ/RawMEA2005/Spikes_all_channels.mat'); toc
% 80 sec from cold; 25 sec the 2nd time.
%spikespy(raw.allElecData(1:10,:));

tic; Ym = mean(raw.allElecData,2);        % subtract mean in lieu of filter
Y = bsxfun(@minus,raw.allElecData, Ym); toc   % 3 sec

wf.d.samplefreq = raw.samplingRate; wf.d.dt = 1/wf.d.samplefreq;

no.meth = 'a';  % since 'j' fails (no empty clips since M>>1)
noi = empiricalnoise(struct('A',Y(:,1:1e5),'samplefreq',wf.d.samplefreq),no)

% empiricalnoise gave eta = 80.8, too big! - why? didn't have zero means!
sopts.verb = 2;
sopts.noi = noi;   % override
[t l] = fit_timeseries(raw.allElecData(:,1:1e4),wf,'Gl',sopts);
% all t-intervals reached 20 rounds of stuffing - bad
% 7952 spikes in 176 s  (40 threads)  for 0.5 sec of data.

% should take 30x for K, 50x for M = 1500x longer per t-pt
% Was for M=7, 1 sec for 8 min

% k=318 found 1918 times!
