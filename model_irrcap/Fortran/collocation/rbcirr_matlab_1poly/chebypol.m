
function  poly=chebypol(x,n,a,b)

% n is order
% [a b] is the range over which the approximation is desired
% x needs to be mapped back into -1 1, the domain of definition
%   for Chebyshev polynomials.

% evaluates the n+1 cheby polynomials of nth order 
% and stacks them in poly= [order 0 at x 
                           %order 1 at x ...]

xrebased = (x -a)/(b-a)*2-1;


poly=zeros(n+1,1);

poly(1) = 1;

if n>0
    poly(2) = xrebased;
end

if n>1
    for i=3:n+1
        poly(i)=2*xrebased*poly(i-1)-poly(i-2);
    end
end

