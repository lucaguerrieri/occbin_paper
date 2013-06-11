% 4th order difference putting NaN as initial values
% so that it does not change the length of vector
% By Matteo Iacoviello

function y = mydiff4(x)

[nrows,ncols] = size(x);

dx=x(5:end) - x(1:end-4) ;

if(nrows > ncols)

    y=[ NaN
        NaN 
        NaN 
        NaN 
        dx];

else

    y=[NaN NaN NaN NaN dx];

end




