function [n x] =get_pdf(z,nbins)


[n x] = hist(z,nbins);
width = x(2)-x(1);
area = sum(n*width);
n = n/area;