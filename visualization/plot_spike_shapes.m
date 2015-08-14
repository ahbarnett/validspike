function plot_spike_shapes(w,title0,vertical_spread,t,fig,o)
% PLOT_SPIKE_SHAPES - plot a set of multi-channel waveforms or clips
%
% plot_spike_shapes(W) where W is a M x Nt x Ns array with M channels, Nt
%  time samples, and Ns spikes (or events) plots a 2D grid of graphs,
%  with the channel increasing vertically downwards, and the spike number
%  increasing horizontally.
%
% plot_spike_shapes(W, title)
% plot_spike_shapes(W, title, vertical_spread)
%  Controls the title for the plot, and the vertical separation between
%  graphs in the signal units (if vertical_spread<0, flip channel ordering to
%  be increasing upwards)
%
% plot_spike_shapes(W, title, vertical_spread, t) adds alignment lines at
%  the times given by t, measured in samples (offset from the first time index).
%  If t is a numeric array, one line per clip is shown.
%  If t is a struct array, the times t(j).t are shown in clip j (ie t is
%   treated as a parameter array p), if present, and/or labels t(j).l are used
%   to label the firing events if present
%
% If W is a struct, it's treated as variable-length clip object with fields
%  X, Ts, tptr, etc (see mergeclips).
%
% plot_spike_shapes(W, title, vertical_spread, t, fig) uses existing figure
%  number.
% 
% plot_spike_shapes(W, title, vertical_spread, t, [], opts) control options:
%    opts.lines = 0,1 switches veritcal dotted lines from t.t
%    opts.nums = 0,1 switches off or on the grey numbers above clips
%
% Without input arguments it runs a self-test using default data

% todo: If the field d is present in W, a 1 ms time bar added
% todo: various tweaks such as allowing axes, figure placement?
% todo: index set input to show subset w/ correct clip # labels.
%
% Magland/Barnett
% 12/16/14 ahb: added self-test, tidied up. 1/23/15 annotation, autosize
% 1/29/15 added optional alignment time plotting. 2/19/15: t may be struct array
% jfm added optional fig parameter (handle to figure to use) - 2/11/15
% 2/19/15 simplified padding mess, and W may be struct for variable-width clips
% 2/27/15 opts, lines, existfig
% 3/31/15 vertical_spread can be <0 for flipping ordering
% 6/10/15 flip sign of vertical_spread, allow labels w/o times


if nargin<1, test_plot_spike_shapes; return; end
if nargin<2 || isempty(title0), title0=''; end
if ~isstruct(w), w = mergeclips(w); end   % convert to struct in variable-len
if nargin<3 || isempty(vertical_spread), vertical_spread = 1.0 * max(abs(w.X(:))); end; % auto scale
if vertical_spread==0, vertical_spread = 1e-16; end % tiny value (in case w.X=0)
vertical_spread = -vertical_spread; % change default to downwards incr
if nargin<4, t=[]; end
if (nargin<5 || isempty(fig)) existfig=0; fh=figure;
else existfig=1; fh=fig; figure(fh); end
if nargin<6, o = []; end               % opts
if ~isfield(o,'lines'), o.lines = 1; end
if ~isfield(o,'nums'), o.nums = 1; end

padding=2;           % means that 1st sample starts at x-coord = padding+1 (=3)
w = padclips(w,padding,nan);
M = size(w.X,1); % spread out the channels
for c=1:M
	if vertical_spread<0, w.X(c,:)=w.X(c,:)+(M-c+1)*abs(vertical_spread);
        else w.X(c,:)=w.X(c,:)+c*vertical_spread; end
end
vertical_spread = abs(vertical_spread);
sx=min(3000,50+50*w.Ns); sy=320;                % controls figure size

plot(w.X'); title(title0);        % show data

if ~existfig, set(fh,'Position', [100,100,sx,sy]); end
xlim([1,w.Ttot]);
ylim(vertical_spread*[.3,M+.8+0.05*M]);

K = w.Ns; % Alex's annotations:
if o.nums && K>1   % number the clips
  for k=1:K, text(w.tptr(k) + w.Ts(k)/2, vertical_spread*(M+0.6), sprintf('%d',k), 'rotation',90,'color',.7*[1 1 1]); end
end
if (length(w.tptr)>1) %condition added by jfm (2/19/15)
	h=vline(w.tptr(2:end), '-'); set(h,'color',.85*[1 1 1]);   % faint dividers
end;
if ~isempty(t)                  % show alignment times
  for k=1:K, p = t(k);          % get kth param struct or scalar
    if ~isstruct(p), tt = p; p = []; p.t = tt; end  % fake a param struct
    if ~isfield(p,'t'), o.lines = 0; p.t = (w.Ts(k)-1)/2*ones(size(p.l));
    end % fake central times
    annotateparams(p, w.tptr(k)+padding, -0.2*vertical_spread, o);
  end
end
set(fh,'Color',[1,1,1]);
set(gca,'YTickLabel',[]); set(gca,'Ytick',[]);  % remove axes
set(gca,'XTickLabel',[]); set(gca,'Xtick',[]);
% qu for jfm: where is box set on?
%%%%

function annotateparams(p,ioff,ytxt,o)
% show parameters p.t (assumed), p.l (if present), offset by index ioff in x.
% ytxt gives text y
if o.lines
  h=vline(ioff + p.t, '--'); set(h,'color',.5*[1 1 1]);
end
if isfield(p,'l')
  for j=1:numel(p.l), h=text(ioff + p.t(j), ytxt, sprintf('%d',p.l(j))); end
end


%%%%%%%
function test_plot_spike_shapes
wf = loaddefaultwaveforms;  % note: are upsampled
W = wf.W;
[M NT K] = size(W); t = ones(1,K) * (NT-1)/2;  % std t_aligns for wf
plot_spike_shapes(W,'test',[],t);
'alignments should be visually correct'
% todo: test variable-length clips
