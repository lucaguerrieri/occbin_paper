% 1st order difference putting 0 as initial values
% so that it does not change the length of vector
% By Matteo Iacoviello

function y = mydiff0(x)

[nrows,ncols] = size(x);

dx=x(2:end) - x(1:end-1) ;

if(nrows > ncols)

    y=[ 0 
        dx];

else

    y=[0 dx];

end




