function writemda(X,fname)
% WRITEMDA - write to JFM's multi-dimensional array file format
%
% writemda(X,fname) writes the array X to the file fname. X can have any number
%  of dimensions up to 6.

% JFM's. todo: finish documentation

num_dims=2;
if (size(X,3)>1) num_dims=3; end;
if (size(X,4)>1) num_dims=4; end;
if (size(X,5)>1) num_dims=5; end;
if (size(X,6)>1) num_dims=6; end;
FF=fopen(fname,'w');
complex=1;
if (isreal(X)) complex=0; end;
if (complex)
    fwrite(FF,-1,'int32');
    fwrite(FF,8,'int32');
    fwrite(FF,num_dims,'int32');
    dimprod=1;
    for dd=1:num_dims
        fwrite(FF,size(X,dd),'int32');
        dimprod=dimprod*size(X,dd);
    end;
    XS=reshape(X,dimprod,1);
    Y=zeros(dimprod*2,1);
    Y(1:2:dimprod*2-1)=real(XS);
    Y(2:2:dimprod*2)=imag(XS);
    fwrite(FF,Y,'float32');
else
    fwrite(FF,-3,'int32');
    fwrite(FF,4,'int32');
    fwrite(FF,num_dims,'int32');
    dimprod=1;
    for dd=1:num_dims
        fwrite(FF,size(X,dd),'int32');
        dimprod=dimprod*size(X,dd);
    end;
    Y=reshape(X,dimprod,1);
    fwrite(FF,Y,'float32');
end;
fclose(FF);
end
