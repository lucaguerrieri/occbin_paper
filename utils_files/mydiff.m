% 1st order difference putting NaN as initial values
% so that it does not change the length of vector
% By Matteo Iacoviello

function y = mydiff(x)

[nrows,ncols] = size(x);

dx=x(2:end) - x(1:end-1) ;

if(nrows > ncols)

    y=[ NaN 
        dx];

else

    y=[NaN dx];

end




