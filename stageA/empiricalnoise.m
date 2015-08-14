function noi = empiricalnoise(d)
% EMPIRICALNOISE - extract noise model parameters from raw EC dataset
%
% noi = empiricalnoise(d). Currently two methods for noi.eta built in.
%
% Inputs:
%  d - raw EC dataset struct, with fields: A - (M*Nt) raw data
%                                          dt - timestep
%                                          etc
% Outputs:
%  noi - noise model struct with fields: eta - std error in iid Gauss model
%
% todo: self-test

% Barnett 2/12/15 

meth='j';   % choose method.   todo: make opt

if strcmp(meth,'a')          % Alex idea: curvature of small signal data
  sc = max(abs(d.A(:))); % fit a Gaussian noise model to small-ampl data...
  b = (-1:0.01:1)'*sc;  % bin centers
  %for m=1:M, %m = 3; h = histc(d.A(m,:),b); % explore each channel separately
  h = histc(d.A(:),b); h = h(:); % all col vecs
  ib = find(h>0.25*max(h)); y = log(h(ib)); x = b(ib); % fit y=const-x^2/(2.eta^2)
  [co] = lscov([1+0*x, x, -x.^2/2],y); eta = 1./sqrt(co(3)); % eta = around 15
  %figure; bar(b,h); set(gca,'yscale','log'); hold on; plot(b,exp(co(1)+co(2)*b-b.^2/(2*eta^2)),'r-');
  noi.eta = eta; % use est std error as Gaussian model

elseif strcmp(meth,'j')      % Jeremy idea: mode of clip-estimated eta distn
  Nc = 1e4;   % # clips
  Tclip = 0.003;  % time clip length in secs. todo: make opt
  T = ceil(Tclip*d.samplefreq); % # samples in time clip
  [M Nt] = size(d.A);
  t = randi(Nt-T,1,Nc); % start indices of clips
  etas = nan(1,Nc);
  for c=1:Nc
    etas(c) = sqrt(sum(sum(d.A(:,t(c):t(c)+T-1).^2)) / (M*T));
  end
  meta = mean(etas); b = meta*(0.5:0.03:10);  % go up to 10x the mean
  f = histc(etas,b);
  %figure; plot(b,f,'+-');
  [~,i] = max(f);      % mode of histogram
  noi.eta = b(i);
end

% todo: meas autocorr & fit to model? This depends on d.A filtering of course
