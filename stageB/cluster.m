function [L cen info] = cluster(z,opts)
% CLUSTER - apply chosen clustering algorithm to feature vectors
%
% L = cluster(z,opts)
% [L cen info] = cluster(z,opts)
%  This is an interface to a variety of clustering methods for pts in R^P.
%  Some algoritm parameters are baked-in here, so can be adjusted.
%  Canonical descending population ordering is used for output clusters.
%
% Inputs
%  z - P x Ns matrix of feature vectors (each column is a point in R^P)
%  opts - controls clustering:
%    opts.cmethod - method:
%         'km' Matlab k-means pre-2014
%         'k' for Jeremy k-means
%         'k++' k-means++
%         'kmed' k-medoids
%         'gmmvb' Gaussian mixture model, variational Bayes, by Mo Chen
%         'dbp' DBSCAN density-based clustering by Patwary (fastest, executable)
%         'dbpe' repeated calls to dbp, auto-adjusting eps
%         'dbd' DBSCAN density-based clustering by Daszykowski
%         'dbk' DBSCAN density-based clustering by Kovesi
%         'frey' Affinity propagation by Frey-Dueck 2007
%         'iso' Magland's isotonic clustering (not implemented; todo)
%    opts.K - number of clusters
%    opts.num_trials - # clustering trials
%                        (w/ different seeds, chooses the best)
%    opts.eps - cluster scale, needed for DBSCAN methods
%    opts.freyp - preference parameter for 'frey', controls K found
%
% Outputs:
%   L - the labels (1 x Ns) - integers between 1 and K
%   cen - the representative cluster centers in feature space (P x K)
%   info - output details of algorithm, including:
%      info.err - the error measure in clustering_err_norm
%      info.fail - success, etc ... todo
%      info.nam - name string of algorithm & # trials
%
% Also see: CLUSTERING_ERR_NRM

% todo: 0) iso 1) hier clus 2) post-2014b mathworks? 3) CURE alg 4) spectral etc
% 5) interface to fastcluster http://danifold.net/fastcluster.html
% 6) interface to pybasicbayes by Matt J Johnson
%
% Barnett 12/19/14. added gmmvb 1/22/15, dbp (and removing 0.7 from eps!) 4/20/15
% 6/12/15: dbpe

if nargin<1, test_cluster; return; end
if nargin<2, opts = []; end
if ~isfield(opts,'K'), opts.K=5; end  % default, arbitrary
if ~isfield(opts,'minpts'), opts.minpts=4; end  % min # pts in clus for DBSCAN
if ~isfield(opts,'cmethod'), opts.cmethod='k++'; end  % default
% Let the particular alg choose the default here...
if isfield(opts,'num_trials'), r=opts.num_trials; else r=[]; end
[P N] = size(z);
K = opts.K;     % requested # clusters
t1 = tic;

% call requested alg (makes L=labels cen=centers), keeping "best" over trials...
if strcmp(opts.cmethod,'km')
  if isempty(r), r=100; end      % default # trials
  info.nam = sprintf('kmeans MathWorks, %d tri (own nrm)',r); % may crash!!
  [L cen] = kmeans_mathworks_preR2014b(z', K,'rep',r);
  cen = cen'; % centroids
  
elseif strcmp(opts.cmethod,'k')
  if isempty(r), r=100; end      % default # trials
  info.nam = sprintf('kmeans Magland, %d tri (own nrm)',r);
  maxit = 100;
  besterr = +inf;
  for i=1:r, [Li ceni erri] = kmeans_jfm(z, K, maxit);
    if erri<besterr, besterr=erri; L = Li; cen = ceni; end   % current best
  end
  
elseif strcmp(opts.cmethod,'k++')
  if isempty(r), r=100; end      % default # trials
  info.nam = sprintf('kmeans++ Sorber, %d tri',r);
  besterr = +inf;
  for i=1:r, [Li ceni] = kmeans_sorber(z, K);
    erri = clustering_err_norm(z,Li,ceni);
    if erri<besterr, besterr=erri; L = Li; cen = ceni; end   % current best
  end

elseif strcmp(opts.cmethod,'kmed')
  if isempty(r), r=100; end      % default # trials
  info.nam = sprintf('k-medoids, Gauss kernel PMTK3, %d tri',r);
  % This is hacked from devel/kmedoidsDemoFaithful.m - ahb:
  noise = 1.0 * opts.eps;      % insensitive over a decade (0.1 fails)
  S = pmtk3_sqdist(z',z'); S = exp(-(1/noise^2)*S); % Gauss kernel
  [med,dpsim] = kmedoids(S,K,r);   % Frey-Dueck... S=similarity matrix
  [score, bestrun] = max(dpsim(end,:)); med = med(:, bestrun); % get best run
  j = unique(med); cen = z(:,j); % indices of medoids in the data, & medoids
  L = zeros(1,N); for i=1:K, L(find(med==j(i)))=i; end % clus labels
  
elseif strcmp(opts.cmethod,'gmmvb')
  if isempty(r), r=10; end      % default # trials
  info.nam = sprintf('Gaussian mixture model, var Bayes, %d tri',r);
  maxK = opts.K;
  bestloglik = -inf;  % use best likelihood (could use err norm)
  for i=1:r
    [Li,modeli,logliks] = vbgm(z, maxK);  % Mo Chen code in GMMvarBayes
    loglik = logliks(end);
    if loglik>bestloglik, bestloglik=loglik; L = Li; cen = modeli.m; end
  end
  Kfound = max(L);  % <=maxK
  
elseif strcmp(opts.cmethod,'dbp')  % deterministic (1 trial ok)
  info.nam = 'dbscan Patwary fast cmdline OpenMP';
  [L Kfound] = dbscanpatwaryexec(z,opts.minpts,opts.eps); % wrapper, disk I/O
  cen = []; % covers K=0 case
  for j=1:Kfound, cen(:,j) = mean(z(:,L==j),2); end % get centroids
  
elseif strcmp(opts.cmethod,'dbpe')  % deterministic (1 trial ok)
  info.nam = 'dbscan Patwary fast cmdline OpenMP, eps-adjust';
  eps = opts.eps;
  for r=1:10
    [L Kfound] = dbscanpatwaryexec(z,opts.minpts,eps); % wrapper, disk I/O
    frac = sum(L>0)/N;
    if frac<0.4, eps = eps*1.2; fprintf('frac class=%.3g, incr eps...\n',frac)
    elseif frac>0.95 || Kfound<=2, eps = eps/1.2; fprintf('frac class=%.3g, K=%d, decr eps...\n',frac,Kfound)
    else, fprintf('eps = %.3g was good in dbpe\n',eps), break; end
  end
  cen = []; % covers K=0 case
  for j=1:Kfound, cen(:,j) = mean(z(:,L==j),2); end % get centroids
  
elseif strcmp(opts.cmethod,'dbd')  % deterministic (1 trial ok)
  info.nam = 'dbscan Daszykowski';  % seems to find less (better) clusters
  % than other two DBSCAN methods. Why? Are params interpreted the same?
  [L type] = dbscan_daszykowski(z', opts.minpts, opts.eps);
  % type says stuff about if in middle or edge
  Kfound = max(L); cen = [];
  for j=1:Kfound, cen(:,j) = mean(z(:,L==j),2); end % get centroids
    
elseif strcmp(opts.cmethod,'dbk')  % deterministic (1 trial ok)
  info.nam = 'dbscan Kovesi';
  [~, L, cen] = dbscan_kovesi(z, opts.eps, opts.minpts);
  Kfound = size(cen,2);
  fprintf('\tKovesi Kfound=%d clusters, %d/%d unclass pts\n',Kfound,...
          numel(find(L==0)),N);

elseif strcmp(opts.cmethod,'frey')  % deterministic (1 trial ok)
  info.nam = 'Frey-Dueck 2007 affinity prop';
  tic; S = -pmtk3_sqdist(z',z');   % dense N^2 elements of similarity matrix
  fprintf('dense S computed in %.g s\n',toc);
  if ~isfield(opts,'freyp'), opts.freyp=10*median(S(:)); end  % p, v big neg hack
  if ~isfield(opts,'sparse'), opts.sparse = 1; end
  if opts.sparse         % sparsify S, to see if apcluster faster; not much...
    nnzperrow = 5;      % how many nonzero elements to keep per row
    for i=1:N, si = S(i,:); [~,j] = sort(si); c = si(j(nnzperrow)); %only nearest
      S(i,:) = si.*(si>=c); end
    S = sparse(S); [i j s] = find(S); S = [i(:),j(:),s(:)]; % pack for apcluster
  end
  fprintf('calling apcluster (p=%.3g)...\n', opts.freyp(1))
  idx = apcluster(S,opts.freyp);  % call their code
  fprintf('apcluster done\n')
  exem = unique(idx); Kfound = numel(exem); L = 0*idx;  % build labels in 1..K
  for k=1:Kfound, L(idx==exem(k)) = k; end
  cen = [];
  for j=1:Kfound, cen(:,j) = mean(z(:,L==j),2); end % get centroids
    
end
fprintf('cluster done in %.3g s\n',toc(t1))

% handling case of variable # clusters found...
if exist('Kfound','var'), K = Kfound; end   % Kfound overrides opts.K, for now

% relabel into canonical population ordering of the K clusters...
clas = L>=1 & L<=K;          % boolean whether each pts classified
L(~clas) = 0; % standardize label for unclassified pts
pops = histc(L,1:K);         % pops of labeled clusters
[~,p] = sort(pops,'descend');
[~,invp] = sort(p); L(clas) = invp(L(clas));  % invert perm to apply to labels
cen = cen(:,p);              % permute centers

if nargout>2, info.err = clustering_err_norm(z,L,cen); end   % if requested
L = L(:)';               % ensure row vec
%%%%%


function test_cluster
% Crude test clustering algs in R^d. Surgery from Alex's clustering sandbox
% 
% todo: perm the cen's & compare against known centroids cs? less important

d = 3; %10;      % dims the pts are in (R^d) - strongly affects DBSCAN eps
eps = 0.07;   % isotropic rand normal scale (SS looks like about 0.05-0.1)
ns = 500*ones(1,5);          % Choose: cluster populations: identical...
%ns = [1500 500 100 20 5];    % ... or widely ranging
nc = numel(ns);     % true # clusters
N = sum(ns); fprintf('test %d clusters (total %d pts) in %d dims...\n',nc,N,d);
rng(1); cs = rand(d,nc);    % true centers: in unit cube, fix seed
x = []; l = [];      % data, and true labels (row vec)
for c=1:nc, x = [x bsxfun(@plus,cs(:,c),eps*randn(d,ns(c)))];
  l = [l c*ones(1,ns(c))]; end
j = randperm(N); x = x(:,j); l = l(j);   % shuffle the pts: l now true labels
plot_labeled_pts(x,l); title('test cluster: true labels'); drawnow;

rng('shuffle'); % restore non-determinism...
% tell it some "helpful" algorithm params (only used by certain algs):
o.K = nc;                 % true # clusters
o.eps = 0.5*eps;  % param for DBSCAN methods using true noise scale

%meths = {'km','k','k++','kmed','gmmvb','dbp','dbd','dbk','frey'}; % all, incl crashable
meths = {'k','k++','kmed','gmmvb','dbp','dbd','dbk','frey'}; % all uncrashable
%meths = {'dbp','dbd','dbk'}; % just DBSCAN methods

for i=1:numel(meths);                      % ----- loop over methods
  o.cmethod = meths{i};
  [L cen info] = cluster(x,o);
  clas = find(L>0);
  p = getbestshuffling(L,l); L(clas) = p(L(clas)); % make best perm of labels
  plot_labeled_pts(x,L); title(info.nam); drawnow;
  [A pk] = labels_similarity(l,L);
  ngood = sum(diag(A)); %numel(find(idx==l))
  fprintf('%-39s: err = %-5.3g  frac good labels = %.3g\n',...
          info.nam,info.err,ngood/N)
  fprintf('\t pops (0...K): '), fprintf('%5d ',histc(L,0:nc)), fprintf('\n')
end                                        % -----
