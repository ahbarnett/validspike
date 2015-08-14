function points_click_waveform(X,z,L,W,nam,t)
% POINTS_CLICK_WAVEFORM - show point cloud with clicking to pop up waveforms
%
% points_click_waveform(X,z,L,W,nam) plots (possibly labeled) points from z,
%  waits for left-clicks which then pop up windows of the signal waveform from
%  X, and waveform from W.  Return (CR) quits it.
%
% points_click_waveform(X,z,L,W,nam,ta) also shows alignment times with X.
%
% Inputs:
%  X = M*T*Ns array of single event windowed data
%  z = Nfea*Ns feature vectors (one col per event)
%  L = 1*Ns integer labels in 1,..,K or 0,..,K (optional, may be empty)
%  W = M*T*K waveforms (optional, may be empty)
%  nam = dataset name (optional)
%  ta = alignment times in original X sample units (optional)
%
% Note: made OBSOLETE by spikespy

% Barnett 1/23/15. fix for unlabeled pts 1/26/15
% moved vline to plot_spike_shapes 1/29/15

if nargin<1, test_points_click_waveform; return; end
[M T Ns] = size(X);
if nargin<3 || isempty(L), L = ones(1,Ns); end
showW = nargin>3 && ~isempty(W);
if nargin<5, nam = ''; end
K = max(L); if size(W,3)<K, warning('not enough waveforms in W for labels'); end

plot_labeled_pts(z,L); f = gcf;
title(sprintf('left-click selects, arrows rotate: %s',nam))
hs = get(gca,'children');
Kh = numel(hs); % # obj in fig
if Kh>K+1, warning('too many objs in fig; I don''t know what they are'); end
% if Kh=K+1 we hope it's due to unlabeled grey + signs

vs = max(abs(W(:))); % vertical spread (don't let it vary by waveform in X)

d = 3;  % # degrees to rotate by using arrow keys
c = 0;   % how many successful clicks
b = 0;   % key press or mouse click
while ~isempty(b)         % "Enter" quits the loop
  b = 0;
  while b~=1    % collect kbd input until get left-click
    [x y b] = ginput(1);     % 2d data from this not actually used
    if b==28, view(get(gca,'view')+[d 0]); drawnow;
    elseif b==29, view(get(gca,'view')+[-d 0]); drawnow;
    elseif b==30, view(get(gca,'view')+[0 -d]); drawnow;
    elseif b==31, view(get(gca,'view')+[0 d]); drawnow;
    end
  end
  k = [];
  if ~isempty(b) && ~isempty(gco)
    k = Kh+1 - find(hs==gco); % cluster #, seems hs comes in reverse order
    if k>K, k=0; end  % unlabeled pt
    [pout, vout, viout, facevout, faceiout]  = select3d(gco); % Conti's code  
  end
  if ~isempty(k) && ~isempty(viout)  % Conti finds something (labeled!)
    i = find(L==k);   % events with current label k (including unlabeled L=0)
    j = i(viout);      % event # selected
    hold on; text(z(1,j),z(2,j),z(3,j),sprintf('%d',j)); % label the z pt
    Y = X(:,:,j);
    if 0&showW && k<=K          % augment w/ avg waveform
      Ytmp = Y; Y = zeros(M,T,2); Y(:,:,1)=Ytmp; Y(:,:,2)=W(:,:,k);
    end
    tit = sprintf('%d k=%d',j,k);
    if nargin<=5, plot_spike_shapes(Y,tit,vs);
    else, plot_spike_shapes(Y,tit,vs,t(j)); end % include alignment times
    fp = get(f,'position');  % move the fig closer to the parent
    g = gcf; gp = get(g,'position'); gp(1:2) = fp(1:2) + [c*100,-400];
    set(g,'Position', gp);
    drawnow
    figure(f);
    c = c+1;
  end
end
%%%%


function test_points_click_waveform   % doesn't test unlabeled pts
fprintf('please wait 20 s while spike sorting done...\n')
d = loaddemodata; o.K = 5; o.quan=0.005;
d.A = freqfilter(d.A,d.samplefreq,300,10000); [M N] = size(d.A);   % Stage A:
[X t m info] = detectevents(d,o); 
o.upsampfac = 4; %1;    % Stage B:
o.realign = 1;
o.fmethod = 'pca'; o.num_fea = 10; % choose one of three feature types
o.cmethod = 'k++';   % choose clustering method
[L W z cinfo ta] = ssalg_feature_cluster(X,o);
plot_spike_shapes(W,'mean waveforms');
points_click_waveform(X,z,L,W,['test ' d.name],ta)
