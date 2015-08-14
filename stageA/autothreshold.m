function t = autothreshold(d)
% AUTOTHRESHOLD  choose threshold for detection given raw EC data struct d
%  containing d.A timeseries, d.samplefreq
%
% function t = autothreshold(d) returns a single threshold based on a heuristic
%  noise model

% Barnett 6/11/15

noi = empiricalnoise(d);
%fprintf('\tnoi.eta = %.3g\n',noi.eta)
t = 5.0*noi.eta;   % how many std deviations from zero
