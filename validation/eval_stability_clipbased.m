function [fhat,fsam,info] = eval_stability_clipbased(alg,X,ssopts,o)
% EVAL_STABILITY_CLIPBASED  stability metrics for clip-based spike sorting alg
%
% [fhat,fsam] = eval_stability_clipbased(alg,X,ssopts,opts)
% [fhat,fsam,info] = eval_stability_clipbased(alg,X,ssopts,opts)
%
% Inputs:
%  alg - function handle to a Stage B clip-based spike sorter with interface:
%   <begin interface>
%   L = alg(Z,ssopts)
%    Inputs to alg:
%     Z is a 3d array (M x Nt x N) of individual spike waveforms
%      M = # of channels              (will be same as in X)
%      Nt = # timepoints per waveform (will be same as in X)
%      N = # spikes                   (may be different from Ns in X)
%     ssopts is as passed to eval_stability_clipbased.
%    Outputs from alg:
%     L = 1-by-N list of labels for clips
%   <end interface>
%
%  X - 3D array (M x Nt x Ns) of event clips to be sorted
%  ssopts - is passed to the SS alg.
%  opts (optional) controls various options:
%   opts.meth - stability metric:
%          'rerun' - vanilla rerunning of SS alg
%              opts.allpairs: if 0 (default), Nr indep pairs of runs,
%                             if 1, Nr runs, uses all Nr(Nr-1)/2 pairs.
%          'cv3' - three-way cross-validation via classifier (2*Nr runs of
%                   size N/3 each)
%          'blur' - self-blurring (Nr+1 runs)
%          'rev' - noise reversal (2 SS runs, ignores Nr)
%
%   opts.num_runs - Nr, controls how many times to run the SS alg
%   opts.verb  - verbosity of output while running (& passed to labels_accuracy)
%
% Outputs:
%    fhat - (1-by-K) best estimates of stability metrics for each neuron type
%    fsam - (Nsam-by-K) complete set of samples of metrics for each neuron type
%    info - struct with:
%           Qs - cell array of confusion matrices produced
%           pops - populations from first full run (rerun,blur,rev only)
%
% See also: DRIVER_CLIPS for example of interface to alg

% simplified from evalconsistency.m. Barnett 6/16/15. blur, noise rev 7/1/15
% CV classifier permute 7/3/15.

% todo: variable-K versions

if nargin<4, o=[]; end
if ~isfield(o,'meth'), o.meth = 'cv3'; end
if ~isfield(o,'verb'), o.verb = 1; end
if isfield(o,'num_runs'), Nr = o.num_runs; else Nr=20; end  % > dozen (Mackay96)
[M Nt N] = size(X);
if o.verb, fprintf('eval_stab_clipbased: # clips = %d (M=%d, Nt=%d)\n',N,M,Nt); end

fsam = []; info.Qs = {};

if strcmp(o.meth,'rerun') % ================== simple rerunning on full set
  if ~isfield(o,'allpairs'), o.allpairs = 0; end
  for r=1:Nr, fprintf('ss run %d\n',r)
    L{r} = alg(X,ssopts);       % do the SS
  end
  K = max(L{1}); pops = histc(L{1},1:K);  % fix K and quoted pops
  if o.verb
    disp('populations in "main" run:'), fprintf('%6d ',pops), fprintf('\n')
  end
  if ~o.allpairs        % plain compare runs in pairs
    for i=1:2:2*Nr-1
      [~,Q,acc] = labels_accuracy(L{i},L{i+1},o);
      fsam = [fsam; acc.p];    % todo *** only works if K fixed
      info.Qs = {info.Qs{:} Q};
    end
  else        % all N(N-1)/2 pairs, won't be statistically indep, but more data
              % todo: handle acc.p different lengths, fix to K for L{1}
    for i=1:Nr
      for j=1:i-1
        [~,Q,acc] = labels_accuracy(L{i},L{j},o);
        fsam = [fsam; acc.p];  % *** only works if K fixed
        info.Qs = {info.Qs{:} Q};
      end
    end
  end

elseif strcmp(o.meth,'cv3') % ============ 3-way CV, with general classifier
  for r=1:Nr, fprintf('ss CV run %d\n',r)
    i = randperm(N); n = round(N/3);
    ia = i(1:n); ib = i(n+1:2*n); ic = i(2*n+1:end);  % index sets
    L = alg(X(:,:,ia),ssopts);       % 1/3 SS
    Ca = build_classifier(X(:,:,ia),L);
    if r==1
      Ca1 = Ca;                              % save 1st classifier, so can...
    else
      Ca = best_permute_classifier(Ca1,Ca);  % ...match later orderings to it
    end
    L = alg(X(:,:,ib),ssopts);       % 1/3 SS
    Cb = build_classifier(X(:,:,ib),L); clear L
    La = classify(X(:,:,ic),Ca);
    Lb = classify(X(:,:,ic),Cb);
    [~,Q,acc] = labels_accuracy(La,Lb,o);
    
    fsam = [fsam; acc.p];  % *** only works if K fixed
    info.Qs = {info.Qs{:} Q};   
  end
  
elseif strcmp(o.meth,'blur') % ============ self-blurring
  if ~isfield(o,'gamma'), o.gamma = 1.0; end
  fprintf('ss blur unperturbed run\n')
  L = alg(X,ssopts);                 % all Nr runs compared to this ordering
  K = max(L);
  W = meanwaveforms(X,L);
  for r=1:Nr, fprintf('ss blur run %d\n',r)
    Xp = X;                                       % will be perturbed
    for k=1:K
      j = find(L==k);
      d = bsxfun(@minus, X(:,:,j), W(:,:,k));     % differences from centroids
      Xp(:,:,j) = Xp(:,:,j) + o.gamma*d(:,:,randperm(numel(j)));    % do blur
    end
    Lp = alg(Xp,ssopts);
    [~,Q,acc] = labels_accuracy(L,Lp,o);
    fsam = [fsam; acc.p];  % *** only works if K fixed
    info.Qs = {info.Qs{:} Q};  
  end
  
elseif strcmp(o.meth,'rev') % ============ noise-reversal
  fprintf('ss noise-rev unperturbed run\n')
  L = alg(X,ssopts);
  K = max(L);
  W = meanwaveforms(X,L);
  Xp = X;                                       % will be noise-flipped
  for k=1:K
    j = find(L==k);
    Xp(:,:,j) = bsxfun(@minus, 2*W(:,:,k), Xp(:,:,j)); % inversion thru centroid
  end
  fprintf('ss noise-rev reversed run\n')
  Lp = alg(Xp,ssopts);
  [~,Q,acc] = labels_accuracy(L,Lp,o);
  fsam = acc.p; info.Qs = Q;  % only one sample!

end                        % =========== done validation algs
fhat = mean(fsam,1);  % best estimate is mean
if iscell(L), L = L{1}; end, info.pops = histc(L,1:K);  % info out

% Report card...
if o.verb
  disp('mean stabilities per spike type:'), fprintf('%6.3f ',fhat), fprintf('\n')
end
if o.verb>1, show_stabilities(fhat,fsam,o); end


%%%%%%%%%%%

function cl = build_classifier(X,L)   % return classifier object, for CV method
cl.W = meanwaveforms(X,L);

function L = classify(X,cl)           % use classifier object, for CV method
[M T N] = size(X);
K = size(cl.W,3);
L = nan(1,N);
for j=1:N
  d = bsxfun(@minus, X(:,:,j), cl.W); % differences from centroids (M*T*K)
  dist = sum(reshape(d,[M*T K]).^2);  % l2 distances (1*K)
  [~,L(j)] = min(dist);
end

function C = best_permute_classifier(C1,C) % permute classifier C to best
% match C1. Particular to the above classifier (which uses mean waveforms).
% Used by CV method. Barnett 7/3/15
% todo: test rectangular case, fix case K1>K
[M T K] = size(C.W); [M T K1] = size(C1.W);
Y = reshape(C.W,[M*T K]); Y1 = reshape(C1.W,[M*T K1]); % sets of real vectors
yy = sum(Y.^2,1); yy1 = sum(Y1.^2,1);  % row vecs of squared norms
D = repmat(yy1',[1 K]) + repmat(yy,[K1 1]) - 2*Y'*Y1;  % squared distances
A = Hungarian(D);
i = zeros(1,K1);            % jfm code from gestbestshuffling_hungarian()
for k=1:K1
  t = find(A(k,:)==1);
  if numel(t)>0, i(k) = t(1); else, i(k) = k; end
end  % only for square case
C.W = C.W(:,:,i);    % permute the C classifier's waveforms
