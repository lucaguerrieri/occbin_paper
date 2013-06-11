RHOA = 0.9;
SIGMA_EPS = 0.01;
%nstates = 81;
%nstates2 = 30;
nstates = 61;
nstates2 = 12;

[P,theta]=markovappr2(RHOA,SIGMA_EPS,3,nstates,nstates2);


[chain1]=markov_match_ar(RHOA,theta(5),theta,50);
[chain2]=markov_match_ar(RHOA,theta(nstates-5+1),theta,50);
chain = [ceil(nstates/2)*ones(1,9) chain1 chain2];

theta_path = theta(chain);
figure
plot(theta_path*100)