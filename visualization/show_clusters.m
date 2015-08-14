function show_clusters(X,L,K)
% SHOW_CLUSTERS - view PCA coords as clusters in 3D plot with colors for the different labels
%
% X - the data (M x Nt x Ns)
% L - the labels (Ns x 1)
% K - the number of clusters
%
% Note that we use pca in this function to compute the first 3x1 feature vectors
%
% Also see: spikespy, which has PCA feature space viewer for clips

% 12/18/14 - ahb changed to use plot_labeled_pts
% 12/19/14 - ahb made use features()

VV=features(X);         % use default features
plot_labeled_pts(VV,L); % Alex plot clusters (shows first 3 dims only)
