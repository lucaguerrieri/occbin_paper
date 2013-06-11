% Transform a transition matrix in its perfect foresight counterpart
% logy: is vector of values for log income
% P: is transition of values for log income
% RHO: is the autoregressive coefficient for log income

function F = perfforappr(P,logy,RHO)

np = size(P,1);
npu = (size(P,1)+1)/2;
npl = (size(P,1)-1)/2;

if abs(sum(sum(P'))-size(P,1))<0.0000001
  disp('Ok, the original matrix is a transition matrix')
  LRP=P^10000;
%   disp(' ')
%   disp('The ergodic probabilities are given by')
%   disp(LRP(npu,:))
else
  error('The original matrix is not a transition matrix')
end


if mod(size(P,1),2)==0
  %number is even
  error('This is a transition matrix with even events, not good for my case')
end

np = size(P,1);
npm = (size(P,1)+1)/2;

F0=0*P;
F1=F0;

for i=1:np
  % For each row
  Elogy=RHO*logy(i);
  iq=max(findnearest(Elogy,logy,0));
  F1(i,iq)=1;
end



if sum(sum(F1'))==size(F1,1)
  disp(' ')
  disp('Ok, the new matrix is a transition matrix')
  LRF=F1^100000;
  disp(' ')
%   disp('The ergodic probabilities are given by')
%   etp=LRF(npm,:);
%   disp(etp)
end

F=F1;

%-----------------------------
% Uncomment below to plot old and new transition matrix
%-----------------------------
% subplot(2,1,1)
% contour(logy,logy,P,100); hold on
% plot(logy,logy,'Linewidth',2)
% title('old transition matrix P')
% subplot(2,1,2)
% contour(logy,logy,F,100); hold on
% plot(logy,logy,'Linewidth',2)
% title('new transition matrix F')



end