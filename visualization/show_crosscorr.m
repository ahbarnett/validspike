function [C taus] = show_crosscorr(l,t,varargin)
% SHOW_CROSSCORR - plot cross-correlation vs time shift given times and labels
%
% [C taus] = show_crosscorr(l,t)
% [C taus] = show_crosscorr(l,t,a)
% [C taus] = show_crosscorr(l,t,a,opts)
%  shows cross-correlation graphs.
%
% Interface is the same as crosscorr.m
% Inputs:
%  l - 1D array of labels
%  t - 1D array of firing times
%  a - 1D array of firing amplitudes (set to 1 if empty or absent)
%  opts controls various options:
%   opts.dtau, taumax = tau time bin width and maximum tau, in sample units.
%   dtau should be odd
% Outputs:
%  C - K*K*Nt matrix of cross-correlations (K = max(l), Nt = # timeshift bins)
%  tau - list of timeshift bins

% Barnett 3/1/15, test and sign-color-info 4/8/15
if nargin<1, test_show_crosscorr; return; end

fprintf('estimating cross-corr...\n')
tic, [C taus] = crosscorr(l,t,varargin{:}); toc
K = size(C,1);
figure;
for i=1:K, for j=1:K, tsubplot(K,K,j+(i-1)*K);
    h = bar(taus,squeeze(C(i,j,:)),'edgecolor','none','barwidth',1.0);
    ch=get(h,'children');     % sign indicated by color
    caxis([-1 1]); set(ch,'FaceVertexCData',-sign(squeeze(C(i,j,:))));
    %axis([-maxtau maxtau sc1 sc2]);
    axis tight off
    %axis tight; set(gca,'xtick',[]);
  end, end
%title('C_{ij}(\tau)'); 
drawnow
%%%

function test_show_crosscorr
pops = [10000 3000 1000 300];
T = 1e6;
K = numel(pops); N = sum(pops);
l = []; for k=1:K; l = [l k*ones(1,pops(k))]; end % make labels
t = T*rand(1,N); % uncorr Poisson spikes: each Cij should be flat
[C taus] = show_crosscorr(l,t);
