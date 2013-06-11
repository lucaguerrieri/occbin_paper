% Forward moving average
% so that it does not change the length of vector
% and for the last periods it only uses limited set of observations
% By Matteo Iacoviello

function y = movavgfwd(x,p)

nx = numel(x);
ma_x = NaN*x;

for i=1:(nx-p+1)

ma_x(i) = mean(x(i:i+p-1)) ;

end

for i=nx-p+1:nx
  ma_x(i) = mean(x(i:end)) ;
end  



y = ma_x ;

% if(nrows > ncols)
% 
%     y=[ NaN 
%         dx];
% 
% else
% 
%     y=[NaN dx];
% 
% end
% 



