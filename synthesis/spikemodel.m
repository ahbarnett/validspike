function Y = spikemodel(wf, p, Nt)
% SPIKEMODEL - forward model for event window given spike labels, times, ampls
%
% Y = spikemodel(wf, p, Nt)
%
% Inputs:
%  wf - waveform object:
%    W     = M*NT*K upsampled centered waveforms (NT assumed odd)
%    fac   = factor to downsample (1, or factor which upsampling was done)
%  p - parameters object:
%    l     = 1*Ns labels in 1,..,K (Ns is the number of spike events)
%    t     = 1*Ns time alignments peak, measured in output samps rel to 1st samp
%    a     = 1*Ns amplitudes (typically close to 1; if empty, all set to 1)
%  Nt    = # output time samples
%
% Outputs:
%  Y     = M*Nt signal window
%
% Obsolete; replaced by: stageC_fitlib/spikemod MEX interface.

% ideas from fitwindow01/spikemodel.m. Barnett 1/28/15

if nargin<1, test_spikemodel; return; end
[M NT K] = size(wf.W);  % NT is # upsampled time points (odd)
if mod(NT,2)==0, warning('waveform NT must be odd'); end
if Nt<1, error('Nt must be positive!'); end
Y = zeros(M,Nt);     % start with zero output before add spikes
if isempty(p), return; end  % no spikes -> zero model signal
Ns = numel(p.l);  % # spikes
if ~isfield(p,'a') || isempty(p.a), p.a = 1+0*p.t; end
if numel(p.t)~=Ns || numel(p.a)~=Ns, error('t and a must be same size as l');end
if wf.fac~=floor(wf.fac), error('wf.fac must be positive integer!'); end

icenW = (NT+1)/2;    % central (peak alignment) index in W grid (NT odd)

for c=1:Ns   % loop over # spikes in param
  iput = 1:Nt;                            % time indices to write to in Y
  iget = round(icenW + wf.fac*(iput - (p.t(c)+1)));  % indices to get from in W
  iok = find(iget>0 & iget<=NT);           % keep only indices inside waveform
  iput = iput(iok); iget = iget(iok);
  Y(:,iput) = Y(:,iput) + p.a(c)*wf.W(:,iget,p.l(c));
end
%%%%

function test_spikemodel
wf = loaddefaultwaveforms;
Nt = 30;

p.t = 14.5; p.l = 1;
Y = spikemodel(wf, p, Nt); plot_spike_shapes(Y,'',[],p.t); % line should hit
p.t = [5.6 14.5]; p.l = [2 3];
Y = spikemodel(wf, p, Nt); plot_spike_shapes(Y,'',[],p.t);

figure; p.l = 1;   % animation sliding ampl and time together
for t=0:0.1:20, p.t = t; p.a = t/20;
  Y = spikemodel(wf, p, Nt); plot(0:Nt-1, Y, '-'); vline(p.t);
  axis([0 Nt-1 -300 150]); drawnow; %pause(0.1);
end
