% figures for time-series validation, EJ 2005 elec359 dataset
% Barnett 7/15/15

clear; d = loaddata('e');  % gets from data_external/
d.A = freqfilter(d.A,d.samplefreq,300,[]);

co.cmethod='k++'; co.K = 9; co.Kkeep = 9; co.thresh=120; co.verb = 2;
so.skip = 5; so.nlps = 10;
[t l p wf R] = spikesort_timeseries(d.A,d.samplefreq,co,[],so);  % do the alg

plot_spike_shapes(wf.W,'waveforms found'); drawnow;
spikespy({d.A,t,l,d.name},{R,'residual'});
show_crosscorr(l,t);  % physiological (refractory hole) check
