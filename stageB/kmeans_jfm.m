function [labels,means,sserr]=kmeans_jfm(vectors,num_clusters,max_iterations)
% KMEANS_JFM - K-means clustering algorithm, implemented by J. Magland.
%
% [labels,means,sserr]=kmeans_jfm(vectors,num_clusters,max_iterations)
%
% Note that sserr is returned but it doesn't use clustering_err_norm()
%
% Pulled from ssalg_pcakmeans_jfm.m  Barnett 12/19/14

	n=size(vectors,1);
	m=size(vectors,2);
	K=num_clusters;
	
	%initialize the cluster centers
	means=kmeans_initialize_cluster_centers(vectors,K);
	
	dists=zeros(m,K);
	for kk=1:K
		means0=repmat(means(:,kk),1,m);
		dists(:,kk)=sqrt(sum((vectors-means0).^2,1))';
	end;
	something_changed=0;
	[minvals,mininds]=min(dists,[],2);
	labels=mininds;
	
	%labels=floor(rand(m,1)*K)+1;
	
	for ii=1:max_iterations
		means=zeros(n,K);
		for kk=1:K
			inds0=find(labels==kk);
			if (length(inds0)>0)
				means(:,kk)=mean(vectors(:,inds0),2);
			end;
		end;
		dists=zeros(m,K);
		for kk=1:K
			means0=repmat(means(:,kk),1,m);
			dists(:,kk)=sqrt(sum((vectors-means0).^2,1))';
		end;
		sserr=0;
		for ii=1:m
			sserr=sserr+dists(ii,labels(ii)).^2;
		end;
		
		something_changed=0;
		[minvals,mininds]=min(dists,[],2);
		inds0=find(labels~=mininds);
		if (length(inds0)>0)
			something_changed=1;
			labels(inds0)=mininds(inds0);
		end;
		if (~something_changed) 
			%disp(sprintf('Ending k-means after %d iterations.',ii));
			return; 
		end;
	end;
	
	disp(sprintf('K-means: Max iterations reached (%d)',max_iterations));
end


function centers=kmeans_initialize_cluster_centers(vectors,num_clusters)

n=size(vectors,1);
m=size(vectors,2);
K=num_clusters;

if (n<K)
	centers=repmat(vectors(:,1),1,K);
	return;
end;

centers=zeros(n,K);
for k=1:K
	centers(:,k)=vectors(:,randi([1,n]));
end;

end
