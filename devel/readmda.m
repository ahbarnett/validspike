function A=readmda(fname)
%READMDA - read the contents of a .mda file. MDA stands for
%multi-dimensional array.
%
% See http://magland.github.io//articles/mda-format/
%
% Syntax: A=readmda(fname)
%
% Inputs:
%    fname - path to the .mda file
%
% Outputs:
%    A - the multi-dimensional array
%
% Other m-files required: none
%
% See also: writemda

% Author: Jeremy Magland
% Jan 2015; Last revision: 15-Feb-2106

F=fopen(fname,'rb');

try
code=fread(F,1,'long');
catch
    error('Problem reading file: %s',fname);
end
if (code>0) 
    num_dims=code;
    code=-1;
else
    fread(F,1,'long');
    num_dims=fread(F,1,'long');    
end;

S=zeros(1,num_dims);
for j=1:num_dims
    S(j)=fread(F,1,'long');
end;
N=prod(S);

A=zeros(S);
if (code==-1)
    M=zeros(1,N*2);
    M(:)=fread(F,N*2,'float');
    A(:)=M(1:2:prod(S)*2)+i*M(2:2:prod(S)*2);
elseif (code==-2)
    A(:)=fread(F,N,'uchar');
elseif (code==-3)
    A(:)=fread(F,N,'float');
elseif (code==-4)
    A(:)=fread(F,N,'short');
elseif (code==-5)
    A(:)=fread(F,N,'int32');
elseif (code==-6)
    A(:)=fread(F,N,'uint16');
elseif (code==-7)
    A(:)=fread(F,N,'double');
end;

fclose(F);