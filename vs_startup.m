function vs_startup
% validspike package MATLAB startup. Barnett 8/13/15
% Call this script from your startup.m script

h = fileparts(mfilename('fullpath'));
addpath(genpath(h))                        % gives access to all subdirs
rmpath(genpath(fullfile(h,'.git')))

%system(sprintf('cd %s; make',h))           % run top-level Makefile

% OPTIONAL:
% set the following to your spikespy location (ie the spikespy/matlab directory):
addpath(genpath(fullfile(h,'../spikespy/matlab')))
