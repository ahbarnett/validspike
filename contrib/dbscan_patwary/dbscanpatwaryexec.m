function [L K] = dbscanpatwaryexec(X, minpts, eps)
% DBSCANPATWARYEXEC - system interface to Patwary's fast DBSCAN OpenMP executable
%
% [L K] = dbscanpatwaryexec(X, minpts, eps)
%
% Inputs:
%  X - d-by-N matrix of real-valued coordinates of N pts in R^d
%  minpts - (int) DBSCAN parameter
%  eps - (float) DBSCAN distance parameter
% Outputs:
%  L - 1-by-N integer list of cluster labels assigned to the points
%       (0 if unassigned)
%  K - (int) number of clusters
%
% NOTE: this is a hack pending building MEX interface

% Barnett 4/20/15

[dims N] = size(X);

tmpin='/tmp/dbscan_patwary_tmp.bin';           % temporary files
tmpout='/tmp/dbscan_patwary_tmp.txt';

[status,txtout]=system('nproc');  % get max # threads (cmdline needs it)
% could use: system('cat /proc/cpuinfo | grep processor | wc -l')
ncpu = str2num(txtout); if isempty(ncpu) || ncpu<1, ncpu=1; end

fid = fopen(tmpin,'w');  % write tmp file in binary Patwary format
params = [N dims];
fwrite(fid,params,'uint32');
fwrite(fid,single(X(:)),'float'); % note loop over pts is fast, over dims slow
fclose(fid);

h = fileparts(mfilename('fullpath'));   % current dir
exec = fullfile(h,'omp_dbscan');  % executable
s = sprintf('%s -i %s -b -m %d -e %.15g -t %d -o %s',exec,tmpin,minpts,eps,ncpu,tmpout);
status = system(s);  % call Patwary exec
if status~=0, error('failed'); end

a = textread(tmpout);  % read a tmp ASCII file in Patwary format
L = a(:,2)';  % row vec
K = max(L);  % number of clusters
