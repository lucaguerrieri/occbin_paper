function [g]=skewness(x)

% find the third central moment

m = mean(x);

n = max(size(x));

m3=0;
for i=1:n
    m3 = m3+(x(i)-m)^3;
end
m3=m3/(n-1);

m2 = var(x);

g = m3/m2^(1.5);