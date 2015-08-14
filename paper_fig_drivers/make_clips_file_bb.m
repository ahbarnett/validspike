% output upsampled+aligned clips X file for figures in validation paper
% Buzsaki data
% Barnett 6/24/15

clear; d = loaddata('bb'); % 3e6 of gp5 CingulateCortex BXRat19
d.A = freqfilter(d.A,d.samplefreq,300,[]);
d = channelprewhiten(d,[],struct('verb',0));  % optional, key for buszaki
noi = empiricalnoise(d)
o.thresh = 120; o.verb = 1;
[X t m info] = detectevents(d,o);
%show_detect(d,t,m,info);
fac = 3; kerpars.Tf = 5;       % upsampling params
X = upsample(X, fac, kerpars);
[X ta] = alignspikes(X, fac);
plot_spike_shapes(X(:,:,1:200),'upsampled clips',-500);
%image_spike_shapes(X);

d.A = []; % don't save timeseries
save data_valid/clips_bb_short_th120_3ms_fac3.mat X fac d
