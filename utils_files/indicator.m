% begin indicator.m
function g = indicator(a, b)
% INDICATOR Returns a function g(x) that determines if x is in [a, b]

% Perform some "one time" error checking on a and b here
g = @(x) mytestfunction(x, a, b);

function y = mytestfunction(x, a, b)
y = (a <= x) & (x <= b);
% end indicator.m