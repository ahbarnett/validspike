function fig = viewraw(d,fig)
% VIEWRAW - plot raw data and electrode picture
%
% viewraw(d) where d is raw data struct opens new figure with signals
% viewraw(d,fig) overplots in figure fig
%
% viewraw with no arguments runs a test.
%
% See also: spikespy package, which has many more features than viewraw

% Barnett 11/14/14. removed addpath 12/17/14

if nargin<1, test_viewraw; return; end

sep = 1.5 * max(abs(d.A(:)));  % vertical separation in signal units
[M N] = size(d.A);
t = (1:N)/d.samplefreq;
j = 1:N; % time indices to show: all of it
if nargin<2, figure; else, figure(fig); clf reset; end
plot(t(j),bsxfun(@plus,d.A(:,j)',sep*(0:M-1))); % signal plot, displaced
xlabel('t (s)'); axis tight;
set(gca,'ytick',(0:M-1)*sep,'yticklabel',num2cellstr(1:M)); % label elec #s
title([d.name sprintf(': M=%d N=%d T=%g',M,N,d.T)]);

if isfield(d,'electrodelocs')
  axes('position',[.85 .85 .15 .15]); % inset
  co = get(gca,'colororder'); Nco = size(co,1); % standard color ordering
  for m=1:M, c = co(mod(m-1,Nco)+1,:); % 1x3 color vector matching signal graphs
    plot(d.electrodelocs(1,m),d.electrodelocs(2,m),'.','color',c); hold on;
    text(d.electrodelocs(1,m),d.electrodelocs(2,m),sprintf('%d',m),'color',c);
  end
  set(gca,'xtick',[],'ytick',[]);
  axis equal tight; v=axis; axis(v + 0.5*[-1 1 -1 1]); % pad view a bit
end

%%%%%
function test_viewraw
d = loaddemodata;
tic, viewraw(d); toc
