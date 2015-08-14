function d=loaddata(s)
% LOADDATA - load certain datasets into an EC struct, supplying various missing info
%
% d = loaddata(s) loads into extracellular struct d the dataset chosen by
%  string s. The output d is in the EC format; see synthesis/loaddemodata.m
%  This is an example function useful in our research. The datasets are not
%  distributed as part of validspike.
%
%  s = 'b': Buzsaki rat gp 5.
%      'bb' : longer version
%      'B' : even longer version, filtered
%      'e': EJ Chichilnisky 2005.
%      'h': Harris d533101 tetrode
%
% loaddata with no arguments defaults to synthesis/loaddemodata
% loaddata with no arguments & no outputs runs self-test

% Barnett 11/14/14. New pointers to data 12/17/14. ss_rootpath 1/26/15
% 6/4/15 no ss_rootpath, and auto-gen default

if nargin<1 && nargout<1, test_loaddata; return; end
if nargin<1, d=loaddemodata; end  % default

h = fileparts(mfilename('fullpath')); % current dir
head = fullfile(h,'../data_external');  % points to external raw data directory
if ~exist(head,'dir')
  error('loaddata: you need to link validspike/data_external to point to some data');
end

% Following are loaders for our research datasets ...write your own loader scripts:

if strcmp(s,'e'), d.name = 'EJ 2005-04-26 elec359';
  load([head '/EJ/2005-04-26_elec359'])  % built-in load of .mat file to workspace
  d.A = data;
  z = [0, exp(2i*pi*(1:6)/6)]; d.electrodelocs = [real(z);imag(z)];
  d.samplefreq = samplingRate;
  
elseif strcmp(s,'b'), d.name = 'Buzsaki CingulateCortex BXRat19 gp5 unfilt';
  fid = fopen([head '/Buzsaki/CingulateCortex_BXRat19_gp5_unfilt_1e6.dat']);
  d.A = fread(fid,[10 inf],'float'); fclose(fid); % read as many rows as exist

elseif strcmp(s,'bb'), d.name = 'Buzsaki CingulateCortex BXRat19 gp5 unfilt';
  fid = fopen([head '/Buzsaki/CingulateCortex_BXRat19_gp5_unfilt_3e6.dat']);
  d.A = fread(fid,[10 inf],'float'); fclose(fid); % read as many rows as exist

elseif strcmp(s,'b1'), d.name = 'Buzsaki CingulateCortex BXRat19 gp1 unfilt';
  fid = fopen([head '/Buzsaki/CingulateCortex_BXRat19_gp1_unfilt_3e6.dat']);
  d.A = fread(fid,[10 inf],'float'); fclose(fid); % read as many rows as exist

elseif strcmp(s,'B')
  load([head '/Buzsaki/CingulateCortex_BXRat19_gp5_filt_9.78e6.mat']);
  % already filtered!
  
elseif strcmp(s,'h'), d.name = 'Harris 2000 tetrode d533101';
  % use FMAToolbox load utils. Extracellular is channels 2-5:
  d.A = LoadBinary([head '/Harris2000/d5331/d533101.dat'], ...
                   'nChannels',8,'channels',2:5)';  % note transpose
  d.samplefreq = 1e4;   % from README.txt or .xml, .nrs
  z = exp(2i*pi*(1:4)/4); d.electrodelocs = [real(z);imag(z)]; % I made this up
  % ...but the 4 signals are so similar it can't matter much.
  
  [M N] = size(d.A); t = (1:N)/d.samplefreq;  % clean: keep only good bits
  j=find((t>26 & t<109.49) | t>110.29);
  d.A = d.A(:,j);
end

if s(1)=='b'  % some parameters common to Buszaki
  M=10; d.samplefreq = 20000;
  % guess locations for 10=site shank from Buzsaki64sp neuronexus:
  x(1:2:M)=(-M/2+1:0)*4-8.5; x(2:2:M-2)=(M/2-2:-1:0)*4+8.5; x(M)=0;
  d.electrodelocs = [x; 20*(M-1:-1:0)]/100; % scale to units of 100um
end

% derived useful stuff general to all data sets...
[d.M N] = size(d.A);
d.dt = 1/d.samplefreq;
d.T = N*d.dt;  % total time

%%%%

function test_loaddata  % load some datasets & view
for s='ebh';
  tic; d = loaddata(s); fprintf('load time = %.3g s\n',toc)
  spikespy({d.A,d.name})
end
