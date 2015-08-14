function n = padclips(m,padding,padval)
% PADCLIPS - pad along the time axis a struct of variable-length clips
%
% B = padclips(A,padding,padval) pads with padval the data A.X in A to give B.X,
%  using "padding" on each end of the time axis of each clip.
%
% See also: plot_spike_shapes, which uses it; mergeclips for struct format

% Barnett 2/19/15
n = m;   % copies over unknown fields, and Ns
n.Ts = m.Ts + 2*padding;
n.Ttot = sum(n.Ts);
n.tptr = [1 1+cumsum(n.Ts(1:n.Ns-1))];
M = size(m.X,1);
n.X = padval*ones(M,n.Ttot);   % start a fresh n.X. JFM's tmp array
for c=1:m.Ns               % could do in one go like subsetclips
  j = 0:m.Ts(c)-1;
  n.X(:,n.tptr(c)+padding+j) = m.X(:,m.tptr(c)+j);  
end
