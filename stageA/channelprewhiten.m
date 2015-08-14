function d = channelprewhiten(d,thresh,o)
% channelprewhiten - spatially pre-whiten a timeseries using C from noise clips
%
% function d = channelprewhiten(d,thresh,o)
%  d is a EC dataset object containing d.A timeseries data, d.samplefreq
%  On output d.A is changed. thresh sets the threshold below which a clip is
%  considered noise; if thresh=[] then a noise level is chosen automatically.
%  o - controls options such as:
%       o.verb - 0,1... verbosity

% todo: self-test
%
% Barnett 6/4/15

if nargin<3, o=[]; end
if ~isfield(o,'verb'), o.verb = 0; end
C = empiricalspacetimecorr(d,thresh,o);
if o.verb, set(gcf,'name','before Q'); end
Q=inv(chol(squeeze(C(:,:,(size(C,3)+1)/2)))');   % matrix to transform so
% rows of d.A are indep and unit-norm
for i=1:size(Q,1), Q(i,:) = Q(i,:)/norm(Q(i,:)); end  % preserve norms of rows
if o.verb, figure; imagesc(Q); axis equal; title('Q channel prewhitening');
  caxis([-1 1]*max(abs(caxis))); colorbar; end
d.A=Q*d.A;     % transform data
if o.verb
  C = empiricalspacetimecorr(d,thresh,o); set(gcf,'name','after Q');
end
