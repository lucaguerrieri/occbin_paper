% Make epsilon

load ..\ks\co2a

YY=AA(2:end)';
XX=[ones(numel(YY),1) AA(1:end-1)' ];
% [B,BINT,R2,RINT,STATS] = regress(YY,XX) ;
B = inv(XX'*XX)*(XX'*YY)
AA0=SIGMABAR^(1-ITEC)*(AA-1);
EPSILON(1)=AA0(1);
RHOA1=B(2);
for t=2:numel(AA0)
    EPSILON1(t)=AA0(t)-RHOA1*AA0(t-1);
end



RHOA2=0.999;
for t=2:numel(AA0)
    EPSILON2(t)=AA0(t)-RHOA2*AA0(t-1);
end


RHOA3=0;
for t=2:numel(AA0)
    EPSILON3(t)=AA0(t)-RHOA3*AA0(t-1);
end

save EPSILON EPSILON1 EPSILON2 EPSILON3