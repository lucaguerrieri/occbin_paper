% generate fake data

% N: states
% T: periods
N=10;
T=50;

de(1)=0;
dq(1)=0;

for i=1:N
  for t=2:50
    x1=1*randn
    x2=0.05*randn
    dq(i,t)=0.5*dq(i,t-1)+x1;
    de(i,t)=0.0042*dq(i,t)^3-0.03*dq(i,t)^2-0.035*dq(i,t)+ 0.046+x2;
  end
end

figure
plot(dq(:),de(:),'.')