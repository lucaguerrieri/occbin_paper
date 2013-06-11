clear

global BETA ALPHA DELTAK GAMMAC PSI PHI iss

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Model Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BETA=0.96;
ALPHA=0.33;
DELTAK=0.10;

RHOA = 0.9;

PSI = 0;        % adjustment cost for capital
SIGMA_EPS = 0.013;

% deterministic steady state
kss = ((1/BETA-1+DELTAK)/ALPHA)^(1/(ALPHA-1));
css = -DELTAK*kss +kss^ALPHA;
iss = DELTAK*kss;

order = 6;  % this is the order of the approximating polynomial
% choose an even number so that the constant for the
% polynomial family matches the steady state


nstates = 51; % this is the number of states in the Markov approximation
% to the AR(1) process for technology

% choose an odd number to include the SS

kstart = .9*kss;
kend = 1.5*kss;


[P,theta]=markovappr(RHOA,SIGMA_EPS,3,nstates);

nkgrid = 100;
kgrid = linspace(kstart,kend,nkgrid) ;
nodes = chebyroots(order+1,kstart,kend);
kgrid = sort([kgrid, nodes']);
nkgrid = length(kgrid);

Z = zeros(nstates,nkgrid);

GAMMAC = 2;
PHI = 0.975;
mat = zeros((order+1)*nstates,1);
load_avec_gammac_2
avec = mat;


for thisk=1:nkgrid
    Z(:,thisk) = rfuncirr(avec,kgrid(thisk),theta,P,kstart,kend);
end

figure

l1=plot((kgrid/kss-1)*100,log10(abs(Z(26,:)/css)),'r--');
set(l1,'LineWidth',2)
hold on

title('Comparing Euler Equation Errors Across Models and Solution Techniques')
ylabel('Absolute Value of Euler Equation Errors, Log 10 Scale')
xlabel('Capital, Percent Deviation from Non-Stochastic Steady State')

xlim([-5 15])   

%%%%%%%%%

% kstart = .7*kss;
% kend = 1.3*kss;
%         
% GAMMAC = 2;
% PHI = 0;
% mat = zeros((order+1)*nstates,1);
% load_avec_gammac_2_phi_0
% avec = mat;
% 
% nkgrid = 100;
% kgrid = linspace(kstart,kend,nkgrid) ;
% nodes = chebyroots(order+1,kstart,kend);
% kgrid = sort([kgrid, nodes']);
% nkgrid = length(kgrid);
% 
% 
% h2=subplot(2,1,2);
% for thisk=1:nkgrid
%     Z(:,thisk) = rfuncirr(avec,kgrid(thisk),theta,P,kstart,kend);
% end
% 
% Zvec = Z(26,:);
% 
% %truncate the Zvec at 1e-14
% Zvec(find(abs(Zvec)<1e-14))=1e-14;
% plot(h2,(kgrid/kss-1)*100,log10(abs(Zvec/css)),'r--','LineWidth',2); hold on
% 
% title('RBC Model without Investment Constraint')
% ylabel('Absolute Value of Euler Equation Errors, Log 10 Scale')
% xlabel('Capital, Percent Deviation from Non-Stochastic Steady State')
% 
% xlim([-5 15])   


currentdir = cd;
cd ../../..

plot_euler_errors_pl

eval(['cd ',currentdir])



