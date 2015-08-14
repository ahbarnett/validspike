% test S (glutton score) computation and locvalidmins MEX interfaces
%
% This is a low-level 1-core test; superceded by fit_timeseries self-test.
% Barnett 5/1/15, switched to default W and std synth funcs 6/10/15

clear; verb = 1;
% synthesis... (could use loaddemodata instead)
wf = loaddefaultwaveforms; [M,T,K] = size(wf.W); d = wf.d;   % setup wf
if verb>1, plot_spike_shapes(W,'W'); 'W_k norms:'
  for k=1:K, w = W(:,1:fac:end,k); norm(w(:)), end, end
N = round(1.0*d.samplefreq);                         % 1 second of time
noi = setup_noisemodel(wf.d,N,25);                   % choose noise std dev eta
firingrate = 100; rates = firingrate*ones(1,K);      % mean rates in Hz
[Y pe] = synth_Poissonspiketrain(wf,N,rates,noi,[],0);

% score S...
fac = 3; tpad = 2; o=[]; o.skip=5; %o.shflags = mod(1:numel(tsh),5)==0; % skip 5
tsh = tpad + (0:floor((N-2*tpad-1)*fac))/fac;    % time shifts for S
tic; S = fillscore(wf,Y,tsh,noi,o);
fprintf('S done in %.3g s = %.3g G pts/sec\n',toc,1e-9*numel(tsh)*M*K*T/fac/toc)
if verb>2, spikespy(S); end % note time-scale here is off by factor 1/fac

% local valid minima of S...
nlps = 10*ones(1,K);   % neg log priors (lambda_l) for each type
tic; [jt l s] = locvalidmins(S,nlps); Nm = numel(jt); % find all neg loc mins
fprintf('%d LVMs (cf %d true spikes) found in %.3g s\n',Nm,numel(pe.l),toc)
p.l = l; p.t = tsh(jt);  % time shifts of LVMs (jt 1-indexed)
if verb, figure; plot(tsh,S,'-'); hold on; plot(p.t,s,'.','markersize',20);
  for i=1:numel(jt), text(p.t(i),s(i)-50,sprintf('%d',l(i))); end
  for i=1:numel(pe.t), text(pe.t(i),0,sprintf('%d',pe.l(i))); end
end
