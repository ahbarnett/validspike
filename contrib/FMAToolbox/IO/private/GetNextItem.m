%GetNextItem - Get next item from a Batch object (iterator mode).
%
% This is a helper class to easily read and parse batch files.
%
%  USAGE
%
%    [b,item] = GetNextItem(b);
%
%    b              batch object

% Copyright (C) 2007 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [b,item] = GetNextItem(b)

% Check number of parameters
if nargin ~= 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help GetNextItem">GetNextItem</a>'' for details).');
end

b.currentItem = b.currentItem+1;
b.currentField = 0;
if b.currentItem > size(b.field,1),
	item = [];
else
	item = b.currentItem;
end
