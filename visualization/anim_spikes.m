function anim_spikes(X,skip)
% ANIM_SPIKES - show all spikes in X array sequentially as animation
%
% anim_spikes(X) animates at most around 100 typical spikes
% anim_spikes(X,skip) animates every one in "skip" spikes

% Barnett 1/6/15
% todo: make interactive viewer as in WaveSorter?
[M Nt Ns] = size(X);
if nargin<2, skip = max(1,round(Ns/100)); end   % default skip
figure;
for j=1:skip:Ns
  v = X(:,:,j)';  % note transpose to get time-fast electrode-slow ordering
  plot(v(:),'.-');
  title(sprintf('spike j = %.-5d of Ns=%d',j,Ns))
  axis([1 M*Nt -1 1])
  drawnow
end
