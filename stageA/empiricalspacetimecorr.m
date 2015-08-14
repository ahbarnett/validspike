function C = empiricalspacetimecorr(d,thresh,o);
% EMPIRICALSPACETIMECORR  estimate "noise" autocorrelation over channels and time
%
% C = empiricalspacetimecorr(d,thresh,o)
%
% Inputs:
%  d - signal data object with M-by-M data d.A, d.dt sample time, etc
%  thresh - (optional) threshold below which counted as noise clip, measured if
%           not given or empty.
%  o - options including:
%      o.verb - if true, make a figure
% Outputs:
%  C - M*M*Nt autocorrelation matrix
%
% dtau is fixed at 1 time sample for now

% Barnett 3/13/15-3/22/15. Changed to function 5/18/15
% todo: make self-test
% todo: speed up - takes 2 sec or so for M=7

if nargin<3, o=[]; end
if nargin<2 || isempty(thresh), thresh = autothreshold(d); end
if ~isfield(o,'Twin'), o.Twin = 0.004; end % window length for below-thresh clips
if ~isfield(o,'verb'), o.verb = 0; end
Nt = ceil(o.Twin/d.dt); if mod(Nt,2)==0, Nt=Nt+1; end % odd, window width
maxtau = (Nt-1)/2; taus = -maxtau:maxtau;  % dtau fixed at 1 sample
[M N] = size(d.A);

fprintf('computing space-time cross-corr... ')
n = 1e4;  % # noise-clip trials
i=1; C = zeros(M,M,Nt);  % cross-corr
while i<=n
  j = randi(N,1);
  if j>maxtau && j<=N-maxtau
    X = d.A(:,j+taus);
    if min(X(:))>-thresh   % it's noise, not an "event"
      for m=1:M
        C(m,:,:) = squeeze(C(m,:,:)) + X(m,1+maxtau)*X; % 1+maxtau is central el
      end
      i = i+1;
      if mod(i,round(n/10))==0, fprintf('%d%% ',round(100*i/n)); end
    end
  end  
end
C = C/n;
fprintf('\n')

if o.verb
  figure; sc1 = min(C(:)); sc2 = max(C(:));
  for i=1:M, for j=1:M, tsubplot(M,M,j+(i-1)*M);
      h = bar(taus*d.dt,squeeze(C(i,j,:)),'edgecolor','none','barwidth',1.0);
      ch=get(h,'children');
      set(ch,'FaceVertexCData',-sign(squeeze(C(i,j,:))));
      axis([-maxtau*d.dt maxtau*d.dt sc1 sc2]); axis off;
      %axis tight; set(gca,'xtick',[]);
    end, end
    title(sprintf('C_{ij}(\\tau) over [-%g,%g] ms',o.Twin*1e3,o.Twin*1e3)); drawnow
end
