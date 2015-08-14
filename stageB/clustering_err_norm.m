function sserr=clustering_err_norm(vectors,labels,means)
% CLUSTERING_ERR_NORM - defines an error norm for a clustering of feature vecs
%
% sserr=clustering_err_norm(vectors,labels,means,K)
%  Inputs:
%   vectors - P x Ns matrix of feature pts (eg z)
%   labels - 1 x Ns list of labels in 1...K
%   means - P x K cluster centers (eg centroids)
%
% from jfm's compute_sserr. Barnett 12/19/14

K = size(means,2);
for k=1:K
  inds=find(labels==k);
  if (length(inds)>0)
    vectors(:,inds)=vectors(:,inds)-repmat(means(:,k),1,length(inds));
  end
end
sserr=sum(vectors(:).^2);
