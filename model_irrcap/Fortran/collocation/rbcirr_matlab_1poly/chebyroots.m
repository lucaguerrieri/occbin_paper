function [nodes]=chebyroots(order,a,b)

% see Judd, Algorithm 6.4 step 1

roots=zeros(order,1);

for k=1:order
    roots(k) = cos((2*k-1)/(2*order)*pi);
end

nodes = (roots+1)*(b-a)/2+a;