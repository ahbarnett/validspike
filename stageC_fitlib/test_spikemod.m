% test spikemode, MEX forward model for spikes
% Taken from synthesis/test_spikemodel
% Barnett 2/13/15
clear; wf = loaddefaultwaveforms;
Nt = 30;

p.t = 14.5; p.l = 1;
Y = spikemod(wf, p, Nt); plot_spike_shapes(Y,'',[],p.t); % line should hit
fprintf('diff from spikemodel = %.3g\n',norm(Y - spikemodel(wf, p, Nt)))
p.t = [5.6 14.5]; p.l = [2 3];
Y = spikemod(wf, p, Nt); plot_spike_shapes(Y,'',[],p.t);
fprintf('diff from spikemodel = %.3g\n',norm(Y - spikemodel(wf, p, Nt)))

figure; p.l = 1;   % animation sliding ampl and time together
for t=0:0.1:20, p.t = t; p.a = t/20;
  Y = spikemod(wf, p, Nt); plot(0:Nt-1, Y, '-'); vline(p.t);
  axis([0 Nt-1 -200 100]); drawnow; %pause(0.1);
end
