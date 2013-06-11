function y = detrendcube(x)
%DETRENDQUAD Remove a cube trend from a vector
%   Y = DETRENDQUAD(X) removes the best straight-line fit linear trend from the
%   data in vector X and returns the residual in vector Y.  If X is a
%   matrix, DETREND removes the trend from each column of the matrix.
%


n = size(x,1);
if n == 1,
  x = x(:);			% If a row, turn into column vector
end
N = size(x,1);

[mmm,mmmm,y] = regress(x,...
                           [ ones(numel(x),1) cumsum(ones(N,1)) cumsum(ones(N,1)).^2 cumsum(ones(N,1)).^3 ]);

