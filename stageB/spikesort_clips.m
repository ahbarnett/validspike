function [L W z cinfo ta finfo] = spikesort_clips(X,opts)
% SPIKESORT_CLIPS - Stage B spike sorter to classify all single event clips
%
% [L W] = spikesort_clips(X,opts) runs a stage-B spike-sorting
%  algorithm using the scheme:
%   upsample -> align -> get features (pts in R^P) -> cluster (classify)
%  The upsampling, alignment, feature choice and clustering algorithms are
%  independent and flexible. By default no upsampling nor alignment is done.
%  The labeling order is chosen in descending norm of the W^k mean waveforms.
%
% [L W z cinfo ta] = spikesort_clips(X,opts) returns more info
%
% Inputs:
%   X - the spike event waveforms (M x Nt x Ns)
%   opts - controls the following:
%     Upsampling & alignment options (default is not upsample nor align):
%       opts.upsampfac - 0 do nothing, 1 aligns no upsamp, >1 aligns upsampled
%                        waveforms using value as upsampling factor
%       opts.kerpars - controls upsampling interpolation kernel, see UPSAMPLE
%       opts.realign - if true, realign to the upsampled W (this is
%                        generated in a 0th pass if not present)
%                        Notes: i) largely obsolete due to stage C;
%                              ii) we align in upsampled data, not the fastest
%       opts.ta      - override alignment times, rather than fit them
%     Feature vector options:
%       opts.fmethod - see FEATURES
%     Clustering options:
%       opts.num_fea - number of feature dimensions to use for clustering
%                      (default 10)
%       All other opts passed to clustering method - see CLUSTER.
%
% Outputs:
%   L - the labels (1 x Ns) - integers between 1 and K
%   W - the representative (possibly upsampled) waveforms (M x Nt x K)
%   z - the feature vectors (num_fea x Ns)
%   cinfo - clustering output info, minimally cinfo.nam w/ its name & reps
%   ta - alignment time offsets of events in sample units (1 x Ns)
%   finfo - feature vector finding into.
%
% Also see for self-test: DRIVER_CLIPS

% Barnett 12/19/14
% upsample & align 12/23/14
% realign 1/12/15, override ta 1/30/15. new name 6/10/15. use meanW 7/1/15

if nargin==0, driver_clips; return; end
if nargin<2, opts = []; end           % (for defaults see feature & cluster)
if ~isfield(opts,'verb'), opts.verb=0; end             % default
if ~isfield(opts,'num_fea'), opts.num_fea=10; end      % "
if ~isfield(opts,'upsampfac'), opts.upsampfac=0; end
if ~isfield(opts,'kerpars'), opts.kerpars = []; end
if ~isfield(opts,'realign'), opts.realign = 0; end
if ~isfield(opts.kerpars,'Tf'), opts.kerpars.Tf=5*(opts.upsampfac>1); end % def

% upsample & align (overwrites X)...
[M T Ns] = size(X);
ta = nan(1,Ns); % dummy peak times
if opts.upsampfac>0
  if opts.realign  % do zero-th pass to get labels, upsampled W for L2-alignment
    o = opts; o = rmfield(o,'realign');   % (made obsolete by fitting)
    [realign.L,realign.W] = ssalg_feature_cluster(X,o);
  else, realign = []; end
  X = upsample(X, opts.upsampfac, opts.kerpars);  % tu not needed
  if isfield(opts,'ta')
    [X ta] = alignspikes(X, opts.upsampfac, realign, opts.ta);
  else, [X ta] = alignspikes(X, opts.upsampfac, realign); end
  ta = ta + opts.kerpars.Tf - 1;  % peak times rel to window starts, in samples
  T = size(X,2);
end

% get features
[z finfo] = features(X,opts);
if size(z,1)<opts.num_fea
  warning('not enough features: lowering num_fea'); opts.num_fea = size(z,1);
end

% cluster in feature space
z = z(1:opts.num_fea,:);        % keep only first requested # feature dimensions
[L cen cinfo] = cluster(z,opts);  % (no use for cen)

W = meanwaveforms(X,L);           % estimated waveforms
K = size(W,3);      % allow clustering alg (not opts) to tell us # clusters

% choose standard label ordering: W norms rather than population sizes
Wnrms = sum(reshape(W,[M*T K]).^2,1);
[~,j] = sort(Wnrms,'descend'); [~,jinv] = sort(j);   % get perm & its inverse
W = W(:,:,j); ii=(L>0 & L<=K); L(ii) = jinv(L(ii));  % apply it

if opts.verb>1, plot_labeled_pts(z,L); title('spikesort clips: z');
  plot_spike_shapes(W); title('spikesort clips: W'); drawnow; end
