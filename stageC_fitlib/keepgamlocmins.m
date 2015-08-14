function [jt l s Nm] = keepgamlocmins(jt,l,s,tsh,gamma)
% KEEPGAMLOCMINS - keep local mins that are gamma-local mins w/ given width gamma
%
% jt - list of indices of times in tsh for incoming and output LVMs
% l - spike types for LVMs
% s - s-values for LVMs (including NLPs)
% tsh - list of time-shifts
% gamma - width in times samples for window [t-gamma,t+gamma]
%
% test with test_sort_filtereddata_synth.m

Nm = numel(jt);
keep = logical(l); % to flag only gamma-loc mins
if 1         % call MEX
  for j=1:Nm, keep(j) = minisgammaloc(jt,l,s,tsh,j,gamma); end
  keep = logical(keep);
else         % do matlab (kept for testing purposes)
for j=1:Nm
  if tsh(jt(j))-tsh(jt(max(1,j-1)))<=gamma || tsh(jt(min(Nm,j+1)))-tsh(jt(j))<=gamma
    ii = find(abs(tsh(jt)-tsh(jt(j)))<gamma);
    keep(j) = s(j)<=min(s(ii));   % test for gam-loc min
end, end
end
jt=jt(keep); l=l(keep); s=s(keep); Nm = sum(keep);
fprintf('keeping %d LVMs\n',Nm)
