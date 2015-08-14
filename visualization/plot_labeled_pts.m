function plot_labeled_pts(z,L,o)
% PLOT_LABELED_PTS - shows points in 2d or 3d space colored by their labels
%
% plot_labeled_pts(z,L) where z is 2-by-N or 3-by-N and L is 1-by-N
%  with entries in 1...K, plots each col of z as a point in R^2 or R^3 as
%  appropriate, with color given by corresponding entry in L. If z has >3 rows,
%  just the first 3 are used.
%  Unclassified pits (L=0) are shows as grey plus signs.
%
% plot_labeled_pts(z,L,opts) controls various options:
%  opts.fig = if true, don't open a new figure (add to current fig)
%  opts.nodotted = if true, don't show dotted lines through origin

% todo: self-test
%
% Barnett 12/18/14. Added grey "+" 12/19/14

if nargin<3, o=[]; end
if ~isfield(o,'fig'), o.fig = 0; end
if ~isfield(o,'nodotted'), o.nodotted = 0; end

warning('off','MATLAB:legend:IgnoringExtraEntries'); %by jfm 2/12/15 -- warning from legnum at bottom

K = max(L);
[dims N] = size(z);
if N~=numel(L), error('length of L must match number of columns of z!'); end
c = 'bgrkmcy';  % colors repeat if >7 labels
if ~o.fig, figure; end

if dims>=3
  for n=1:K, j=L==n; % for some obscure reason an empty plot3 locks to 2d view!
    if sum(j), plot3(z(1,j),z(2,j),z(3,j),[c(mod(n-1,numel(c))+1) '.'], 'markersize',5);
      hold on; end % show labels via color
  end
  j=L==0; plot3(z(1,j),z(2,j),z(3,j), '+', 'color',[.5 .5 .5]);
  axis equal vis3d; xlabel('z_1');ylabel('z_2');zlabel('z_3');

elseif dims==2
  for n=1:K, j=L==n; plot(z(1,j),z(2,j),[c(mod(n-1,numel(c))+1) '.']);
    hold on; end % show labels via color
    j=L==0; plot(z(1,j),z(2,j), '+', 'color',[.5 .5 .5]);
  axis equal; xlabel('z_1');ylabel('z_2');
  
else, error('size(z,1) must be 2 or more!');
end

if ~o.nodotted, hline(0); vline(0); end % makes origin obvious
legnum(1:K);     % label colors
