% test for spikefitlib/fitonesp MEX. Tweaked from fitonespike.m
% Barnett 2/13/15. 4/30/15 testing locflag. defaultwaveforms 6/8/15
clear; wf = loaddefaultwaveforms; d = wf.d;
wf.W = wf.W/max(abs(wf.W(:)));  % L-infty normalize
[M,T,K] = size(wf.W);
Nt = 30;
pe.t = 13.35; pe.l = 2;       % pick time and identity
Y = spikemodel(wf, pe, Nt);
noi = setup_noisemodel(d,Nt,0.1);   % must be iid for MEX
o.locflag = 1;      % 0 for regular NLL, 1 for faster one in long clips
[p Jbest f] = fitonesp(wf,Y,noi,o);        % MEX
fprintf('no noise:\n'); pe, p
[p1 Jbest1 f1] = fitonespike(wf,Y,noi);  % matlab
fprintf('differences from matlab: %.3g, %.3g\n',abs(Jbest-Jbest1),norm(f-f1))
