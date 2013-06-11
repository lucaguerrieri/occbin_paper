function [Tran,s]=markovappr2(lambda,sigma,m,N,N2)
% function [Tran,s,probst,alambda,asigmay]=markovappr(lambda,sigma,m,N)
%
% the simple case of approximating first-order 
% autoregressive process with Markov chain
%
% y_t = lambda * y_(t-1) + u_t
%
% u_t is a Gaussian white noise process with standard deviation sigma.
%
% m determines the width of discretized state space, Tauchen uses m=3
% ymax=m*vary,ymin=-m*vary, ymax and ymin are two boundary points
%
% N is the number of possible states chosen to approximate
% the y_t process, usually N=9 should be fine
%
% N2 is the subset of nodes placed between the interval [-sigma sigma]
%
% Tran is the transition matrix of the Markov chain
%
% s is the discretized state space of y_t
%
% alambda is the theoretical first order autoregression coefficient 
% for Markov chain
%
% asigma is the theoretical standard deviation for Markov chain Y_t

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% discretize the state space

stvy = sqrt(sigma^2/(1-lambda^2)); % standard deviation of y_t

ymax = m*stvy   ;                  % upper boundary of state space
ymin = -ymax     ;                 % lower boundary of state space

w = (ymax-ymin)/(N-N2-1);             % length of interval 

s1 = ymin:w:ymax;                   % the discretized state space

ymin2_pos = max(find(s1<-2*sigma));
ymax2_pos = min(find(s1>2*sigma));

s2 = linspace(s1(ymin2_pos),s1(ymax2_pos),N2+2);

s = sort([s1,s2(2:end-1)]);



% calculate the transition matrix

for j=1:N;
   
   for k=2:N-1;
      
      Tran(j,k)= normcdf(s(k)-lambda*s(j)+(s(k+1)-s(k))/2,0,sigma)...
         - normcdf(s(k)-lambda*s(j)-(s(k)-s(k-1))/2,0,sigma);
      
   end
   
   Tran(j,1) = normcdf(s(1)-lambda*s(j)+w/2,0,sigma);      % only subtract half the interval on the right
   Tran(j,N) = 1 - normcdf(s(N)-lambda*s(j)-w/2,0,sigma);  % only subtract half the interval on the left
   
end


   
   
   






