function p=normcdf(x,mu,sigma)

p=(1+erf((x-mu)/sqrt(2*sigma^2)))/2;