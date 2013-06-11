BETA=0.96;
ALPHA=0.33;
DELTAK=0.10;
RHOA = 0.9;
PSI = 0;        
SIGMA_EPS = 0.013;

% steady states
kss = ((1/BETA-1+DELTAK)/ALPHA)^(1/(ALPHA-1));

kap = linspace(lbound*kss,ubound*kss,npoints)';
klog  = -0.05:0.01:0.15;
kgrid = exp(klog+log(kss));
nkgrid = length(kgrid);
    
welf_dp = zeros(1,nkgrid);
for kindx =1:nkgrid
        
        k0 = kgrid(kindx);
        k0pos = max(find(min(abs(kap-k0))==abs(kap-k0)));
    
        kgrid_dp(kindx) = kap(k0pos)   
end
    
    
  