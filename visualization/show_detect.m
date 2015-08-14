function show_detect(d,t,m,info)
% SHOW_DETECT - overlay detection windows (single, multiple) on signal plot
%
% show_detect(d,t,m,info)
% d - raw EC data object
% t - list of window start times as in detectevents
% m - struct containing multiple-spike events
% info - as output by detectevents
%
% Shows single-tagged events w/ pink rectangle, multiple-ones w/ green.

% todo: self-test? not worth it
% Barnett 2/10/15. taken from test_detectevents.

Nt = info.Nt;
Ns = numel(t);  % # single events
viewraw(d);    % time axis in physical units
ax = get(gcf,'children'); set(gcf,'currentaxes',ax(2)); % goto other axes
hold on; v = axis; ylo = v(3); yhi = v(4); % get y range for vertical bars
patch([t;t+Nt;t+Nt;t]*d.dt, repmat([yhi;yhi;ylo;ylo],[1 Ns]), ...
      [1 .8 .8], 'linestyle','none','facealpha',.3);    % single-spike

for i=1:numel(m.t), t = [m.t(i) m.tlast(i)];          % multi-spike [start stop]
  patch([t(1);t(2);t(2);t(1)]*d.dt, [yhi;yhi;ylo;ylo], ...
      [.8 1 .8], 'linestyle','none','facealpha',.3);
end

sc = 0.1*(yhi-ylo)/info.thresh; % also show threshhold criterion
N = size(d.A,2);
hold on; plot((1:N)/d.samplefreq, info.sig*sc,'k-'); hline(info.thresh*sc)
