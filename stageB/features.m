function [z finfo] = features(X,o)
% FEATURES - return feature vectors of chosen type from single events 3d array
%
% z = features(X,opts) where X is M x Nt x Ns 3d single-event data array returns
%  z, a K x Ns array where each column is an event feature vector in R^K, where
%  K is set by the feature type. opts
%  controls the type of feature extracted, thus:
%   opts.fmethod =
%        'pca' PCA applied to signal vector across channels and times (default)
%        'pcachan' PCA applied separately to each channel (across times)
%        'cen' centroid weighting (Prentice et al 2011)
%        'raw' simply the raw data with one column per event (Id operation)
%   opts.d : EC data struct d for electrode info (not needed: d.A, is huge) 
%
% [z finfo] = features(X,opts) also returns finfo.subspace if 'pca', giving
%  a M x Nt x Nfea array of waveform spanning set for PCA projection.
%
% Notes: 1) the returned number of features K is as large as possible; the
%  user can project into a lower-dimensional space later.

% todo: make a data test not dependent on data_external
%
% Barnett 12/19/14. 'cen' 12/22/14, finfo 1/26/15, 'pcachan' 4/17/15

if nargin<1, test_features; return; end
if nargin<2, o = []; end
if ~isfield(o,'fmethod'), o.fmethod='pca'; end  % default
[M T Ns] = size(X);
tic; finfo = [];
% call requested alg...
if strcmp(o.fmethod,'pca'), [z U] = features_pca(X);
  finfo.subspace = reshape(U,[size(X,1),size(X,2),size(U,2)]);
elseif strcmp(o.fmethod,'pcachan')
  z = zeros(M*T,Ns); U = zeros(M*T,M*T);
  for m=1:M                % PCA for each channel separately
    r = m+(0:T-1)*M;  % row indices to write to: all 1st dims together, etc.
    [z(r,:) U(r,r)]= features_pca(X(m,:,:));
  end
  finfo.subspace = reshape(U,[size(X,1),size(X,2),size(U,2)]);
elseif strcmp(o.fmethod,'cen'), z = features_spatial_centroid(X,o.d);
elseif strcmp(o.fmethod,'raw')
  z = reshape(permute(X,[2 1 3]),[M*T Ns]); % matrix: events are cols
else, error('unknown fmethod in features!');
end
fprintf('features done in %.3g s\n',toc)
%%%%%


function [z U] = features_pca(X)
% FEATURES_PCA - get principal component analysis feature vectors, for Ns large
% Jeremy's version using X X^T, faster for large Ns, but limits to 8 digits?
% Assumes that M*Nt < Ns otherwise it's slower than plain SVD on X.
% Tries to standardize signs of z.
% todo: investigate QR or LQ for highly fat matrix SVD case.

[M Nt Ns] = size(X);           % Get some dimensions
MM=M*Nt; X=reshape(X,MM,Ns);   % collapse channel and time dimensions
[U,D] = eig(X*X');   % takes O(MM^2*(MM+Ns)). Note eig faster than svd.
[d,I] = sort(diag(D),'descend'); U = U(:,I);  % sort eigenvectors
U = bsxfun(@times, U, std_signs(U));   % std signs of col vecs
z = U'*X;   % get all components in O(MM^2*Ns).   sing vals = sqrt(diag(D))
% sqrt(d(1:10)), U(1:10,1)   % few singular values & 1st left vec

function z = features_pca_crude(X)
% FEATURES_PCA_CRUDE - get principal component analysis feature vectors.
% Alex's version, plain SVD on A, slower for large Ns, full e_mach.
% Tries to standardize signs of z.
[M T Ns] = size(X);
X = reshape(permute(X,[2 1 3]),[M*T Ns]); % matrix: events are cols
% (the permute has no effect for PCA but is my standard way to unfold to matrix)
[U S V] = svd(X,'econ');     % O((M*T)^2 * Ns) ? not sure for a fat matrix
S = bsxfun(@times, S, std_signs(U));   % std signs of col vecs
z = S*V';  % silly, could use repmat to make O(M*T*Ns)

function s = std_signs(U)
% standardized signs from col vecs of U. Barnett 12/23/14
[m n] = size(U);
s = nan(1,n);
for j=1:n
  [~,i] = max(abs(U(:,j)));   % index of max entry in jth col vec of U
  s(j) = sign(U(i,j));
end

function z = features_spatial_centroid(X,d)
% electrode location features, idea from Prentice et al 2011. Barnett 12/22/14
% Inputs:
%  X - usual 3D single-event array
%  d - raw EC data struct, for electrode loc info
% Outputs:
%  z - feature vectors as 3*Ns array
[M T Ns] = size(X);
z = nan(3,Ns);  % allocate output
for i=1:Ns         % loop over events...
  w = max(-X(:,:,i),[],2);  % weight by negative voltage peak - Prentice
  z(3,i) = max(w);          % 3rd coord is max peak height
  w = w(:)/sum(w);        % normalized weights col vec
  z(1:2,i) = d.electrodelocs * w;  % weighted spatial electrode locs in xy
end


%%%%%%%%%%
function test_features % some basic tests on feature extraction methods
% Barnett 12/19/14. std signs & real data added 12/22/14

if 1 % PCA features for random data
Ns = 1e4; X = randn(3,5,Ns); % Ns large
z = features_pca(X); z2 = features_pca_crude(X); % compare 2 PCA methods
disp('Alex vs Jeremy: differences in each PCA component (signs std)...')
for i=1:size(z,1), disp(norm(z(i,:)-z2(i,:))), end
%figure; imagesc(z); figure; imagesc(z2);
end

if 0 & exist('data_external','dir')  % try features of real datasets
  d = loaddata('b'); disp(d.name)
  d.A = freqfilter(d.A,d.samplefreq,300,[]);
  o.quan = .01; [X t m info] = detectevents(d,o); % old quantile-based detection
  fac = 1 % for 'e' makes little difference
  if fac>1, kerpars.Tf = 5; [X tu] = upsample(X,fac,kerpars); end % upsample
  [Y ta] = alignspikes(X,fac); [M T Ns] = size(Y);
  o.fmethod = 'pca';                      % fea meth: PCA
  z = features(Y,o);
  figure; plot3(z(1,:),z(2,:),z(3,:),'.'); axis vis3d equal
  xlabel('z_1'); ylabel('z_2'); zlabel('z_3'); title('features: PCA')
  o.fmethod = 'cen'; clear d.A; o.d = d;  % fea meth: Centroid spatial
  z = features(Y,o);
  figure; plot3(z(1,:),z(2,:),z(3,:),'.'); axis vis3d %equal
  xlabel('z_1'); ylabel('z_2'); zlabel('z_3'); title('features: elec centroid')
  hold on; plot(d.electrodelocs(1,:),d.electrodelocs(2,:),'r+');
end
% conclusion: centroid features not as clean as PCA ones
% idea: normalize & combine them together to one vector?
