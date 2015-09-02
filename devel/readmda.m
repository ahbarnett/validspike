function A=readmda(fname)

F=fopen(fname,'rb');

code=fread(F,1,'long');
if (code>0) %this is the rare (old) case of complex numbers
    num_dims=code;
    code=-1;
else
    fread(F,1,'long'); % number of bytes per data entry, we skip it
    num_dims=fread(F,1,'long');    
end;

S=zeros(1,num_dims);
for j=1:num_dims
    S(j)=fread(F,1,'long'); %Size of each dimension
end;
N=prod(S);

A=zeros(S);
if (code==-1)
    M=zeros(1,N*2);
    M(:)=fread(F,N*2,'float32');
    A(:)=M(1:2:prod(S)*2)+i*M(2:2:prod(S)*2); %Split into real/imag
elseif (code==-2)
    A(:)=fread(F,N,'uchar');
elseif (code==-3)
    A(:)=fread(F,N,'float32');
elseif (code==-4)
    A(:)=fread(F,N,'int16');
elseif (code==-5)
    A(:)=fread(F,N,'int32');
elseif (code==-6)
    A(:)=fread(F,N,'uint16');
end;

fclose(F);
