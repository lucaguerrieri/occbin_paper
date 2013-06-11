% solves the RBC model by orthogonal collocation with
% Chebyshev polynomials


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


% discretize the asset space


% order = 6;  % this is the order of the approximating polynomial
%             % choose an even number so that the constant for the
%             % polynomial family matches the steady state
%
% % get markov approximation to AR process
% nstates = 61; % this is the number of states in the Markov approximation
%               % to the AR(1) process for technology
%
%               % choose an odd number to include the SS
%
% nstates2 = 12;
  generate_data = 1;
  option = 1;    % set to 1 for random sample
                 % set to 2 for up and down shocks
                 % set to 3 for sequence of technology at steady state 
                 % set to 4 for up shock 
  
  compute_euler_residuals = 1;
  
  compute_welfare=0;
  euler_error_mesh = 0;
  compute_decision_rules = 0;
  load_from_fortran_program = 1;
for version = 1:1
    
    if version ==1
        
        
        % also choose k =3.44 z=-6%
        % and choos k +10% z = 4%
        order = 6;  % this is the order of the approximating polynomial
        % choose an even number so that the constant for the
        % polynomial family matches the steady state
        
        
        nstates = 51; % this is the number of states in the Markov approximation
        % to the AR(1) process for technology
        
        % choose an odd number to include the SS
        
        kstart = .9*kss;
        kend = 1.5*kss;
        
        GAMMAC = 2;
        PHI = 0.975;
        
        file_name_welf = 'welf_results_collocation_gammac_2_ss';
        
    elseif version==2
        order = 6;  % this is the order of the approximating polynomial
        % choose an even number so that the constant for the
        % polynomial family matches the steady state
        
        
        nstates = 51; % this is the number of states in the Markov approximation
        % to the AR(1) process for technology
        
        % choose an odd number to include the SS
        
        kstart = .9*kss;
        kend = 1.5*kss;
        
        GAMMAC = 3;
        PHI = 0.975;
        
        file_name_welf = 'welf_results_collocation_gammac_3_ss';
    elseif version == 3
        order = 6;  % this is the order of the approximating polynomial
        % choose an even number so that the constant for the
        % polynomial family matches the steady state
        
        
        nstates = 51; % this is the number of states in the Markov approximation
        % to the AR(1) process for technology
        
        % choose an odd number to include the SS
        
        kstart = .9*kss;
        kend = 1.5*kss;
        
        GAMMAC = 4;
        PHI = 0.975;
        
        file_name_welf = 'welf_results_collocation_gammac_4_ss';
        
    elseif version == 4
        order = 6;  % this is the order of the approximating polynomial
        % choose an even number so that the constant for the
        % polynomial family matches the steady state
        
        
        nstates = 51; % this is the number of states in the Markov approximation
        % to the AR(1) process for technology
        
        % choose an odd number to include the SS
        
        kstart = .9*kss;
        kend = 1.5*kss;
        
        GAMMAC = 5;
        PHI = 0.975;
        
        file_name_welf = 'welf_results_collocation_gammac_5_ss';
        
        elseif version == 5
        order = 6;  % this is the order of the approximating polynomial
        % choose an even number so that the constant for the
        % polynomial family matches the steady state
        
        
        nstates = 51; % this is the number of states in the Markov approximation
        % to the AR(1) process for technology
        
        % choose an odd number to include the SS
        
        kstart = .7*kss;
        kend = 1.3*kss;
        
        GAMMAC = 2;
        PHI = 0;
        
        file_name_welf = 'welf_results_collocation_gammac_1_phi_0_ss';
        
        elseif version == 6
        order = 6;  % this is the order of the approximating polynomial
        % choose an even number so that the constant for the
        % polynomial family matches the steady state
        
        
        nstates = 3; % this is the number of states in the Markov approximation
        % to the AR(1) process for technology
        
        % choose an odd number to include the SS
        
        kstart = .9*kss;
        kend = 1.5*kss;
        
        GAMMAC = 10;
        PHI = 0.975;
        SIGMA_EPS = 0.03;
        
        file_name_welf = 'welf_results_collocation_gammac_10';
    elseif version == 7
        order = 6;  % this is the order of the approximating polynomial
        % choose an even number so that the constant for the
        % polynomial family matches the steady state
        
        
        nstates = 51; % this is the number of states in the Markov approximation
        % to the AR(1) process for technology
        
        % choose an odd number to include the SS
        
        kstart = .9*kss;
        kend = 1.5*kss;
        
        GAMMAC = 1.01;
        PHI = 0.975;
        
    end
    
    nstates2 = 0;
    
    nnegs = nstates;
    
 
    
    if nstates>1
        %[theta, P] = rouwenhorst(RHOA,SIGMA_EPS,nstates);
        
        [P,theta]=markovappr(RHOA,SIGMA_EPS,3,nstates);
        
        
    else
        P = 1; theta=0;
    end
    
    if ~load_from_fortran_program
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  SOLVE MODEL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % first guess -- keep K constant at its SS value, regardless of theta
        % the zeroth order coefficient is fixed at the SS value
        % fix the first order coefficient at a small negative value -- if capital
        % is above SS, capital is expected to be lower next period.
        %
        % first solve the deterministic model without capital irreversibility,
        % then apply the solution from the deterministic model to the stochastic model
        % without capital irreversibility
        % and finally solve the model with capital irreversibility.
        
        
        avec = zeros(order+1,1);  % first guess for Chebyshev coefficients
        
        % first guesses for the deterministic model
        avec(1) = css^(-GAMMAC);
        avec(2) = -0.1;
        
        % set Markov prob matrix and state vector for deterministic case
        Pd = 1; thetad = 0;  % one state of nature, set at steady state for technology
        
        nodes = chebyroots(order+1,kstart,kend);
        
        rfuncirr(avec,nodes,thetad,Pd,kstart,kend)
        
        % options for fsolve
        options = optimset('Display','Iter','MaxFunEvals',1e10,'MaxIter',1e5,'TolFun',1e-7,'Algorithm','trust-region-reflective');
        
        
        % test residual function
        % rfunc(avec,nodes,theta,P,kstart,kend)
        
        % find deterministic solution of model without capital irreversibility constraint
        A_sol_deterministic = ...
            fsolve(@(avec) rfuncirr(avec,nodes,thetad,Pd,kstart,kend), avec,options);
        
        nnegs = length(find(theta<0));  % number of states in which technology is below SS
        % assume and then verify that
        % capital constraint will not be binding
        % when tehcnology is above SS
        
        
        avec = repmat(A_sol_deterministic,nstates,1); % stacks nstates copy of the deterministic solution
        % find stochastic solution
        A_sol = fsolve(@(avec) rfuncirr(avec,nodes,theta,P,kstart,kend), avec,options);
        
        
        
    else
        
        mat = zeros((order+1)*nstates,1);
        if version == 1
            load_avec_gammac_2
        elseif version == 2
            load_avec_gammac_3
        elseif version ==3
            load_avec_gammac_4
        elseif version == 4
            load_avec_gammac_5
        elseif version == 5 
            load_avec_gammac_2_phi_0
        elseif version == 6
            load_avec_test
        end
        A_sol = mat;
        
    end
    
    avec = A_sol;
    A_sol = reshape(A_sol(1:(order+1)*nstates),order+1,nstates);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot model variables for a given markov chain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
  
    if generate_data
        
        
        % generate path according to markov chain defined below
        %chain_pos1 = markov_match_ar(RHOA,2*SIGMA_EPS,theta,1);
        %chain_pos2 = markov_match_ar(RHOA,-2*SIGMA_EPS,theta,1);
        %chain = [ones(9,1)*ceil(nstates/2) ; chain_pos1*ones(10,1); chain_pos2*ones(10,1); ceil(nstates/2)*ones(10,1)];
        
        
         if option ==1  
         s=RandStream('mt19937ar','seed',1);
         RandStream.setDefaultStream(s);
         nperiods = 200;
         chain=markov_luca(P,nperiods,floor(nstates/2));  
        
         elseif option == 2
             
         [chain1]=markov_match_ar(RHOA,theta(5),theta,50);
         [chain2]=markov_match_ar(RHOA,theta(nstates-5+1),theta,50);
         chain = [ceil(nstates/2)*ones(1,9) chain1 chain2];
         elseif option == 3
         chain = ceil(nstates/2)*ones(1,100);
         else
         [chain1]=markov_match_ar(RHOA,theta(5),theta,50);
         chain = [chain1];
         end
        
       
        nperiods = length(chain);
        k_path = zeros(nperiods+1,1);
        c_path = zeros(nperiods+1,1);
        lambdak_path = zeros(nperiods+1,1);
        resid_path = zeros(nperiods+1,1);
        
        k_path(1) = kss;
        lambdak_path(1) = 0;
        for period = 2:nperiods+1;
            waitbar(period/(nperiods+1))
            state = chain(period-1);
            a_state = [A_sol(:,state)]';
            efunc = a_state*chebypol(k_path(period-1),order,kstart,kend);
            
            
            lambdak_path(period-1) = 0;
            
            c_path(period-1) = efunc^(-1/GAMMAC);
            k_path(period) = invcfunc(theta(state),k_path(period-1),c_path(period-1));
            
            % next verify guess for lambdak
            if k_path(period) -(1-DELTAK)*k_path(period-1)<PHI*iss
                k_path(period) = (1-DELTAK)*k_path(period-1)+PHI*iss;
                c_path(period-1)=cfunc(theta(state),k_path(period-1),k_path(period));
                lambdak_path(period-1) = c_path(period-1)^(-GAMMAC)-efunc;
            end
            
            if compute_euler_residuals
                resids = rfuncirr(avec,k_path(period-1),theta,P,kstart,kend);
                resid_path(period-1) = resids(chain(period-1))/c_path(period-1);
            end
        end
        
        theta_path = theta(chain)';
        
        %%
        figure
        subplot(3,2,1)
        plot(real((k_path(2:end)-kss))/kss*100)
        ylabel('Percent Dev. from SS')
        title('Capital')
        
        subplot(3,2,2)
        i_path = k_path(2:end)-(1-DELTAK)*k_path(1:end-1);
        plot(100*real((i_path-iss))/iss)
        ylabel('Percent Dev. from SS')
        title('Investment')
        
        subplot(3,2,3)
        plot(lambdak_path(1:end-1))
        ylabel('Level')
        title('Lagrangian Multiplier on Kuhn Tucker Condition')
        
        subplot(3,2,4)
        plot((c_path(1:end-1)-css)/css*100)
        ylabel('Percent Dev. from SS')
        title('Consumption')
        
        subplot(3,2,5)
        plot(theta_path*100)
        ylabel('Percent Dev. from SS')
        title('Technology')
        
        if compute_euler_residuals
            subplot(3,2,6)
            plot((resid_path(1:end-1))*100)
            title('Euler Residual (percent of SS consumption)')
        end
        
        % bridge between Luca and Matteo's notation
        
        simC = c_path(1:end-1)';
        simK = k_path(2:end)';
        simI = i_path';
        simZ = exp(theta_path(1:end)');
        simLAMBDAK = lambdak_path(1:end-1)';
        resid_path = resid_path(1:end-1)';
        
        save ../../../luca_nl_sim simC simK simI simZ simLAMBDAK resid_path
    end
    
    
    
    
    if euler_error_mesh
        
       nkgrid = 100;
       Z = zeros(nstates,nkgrid);
      
       kgrid = linspace(kstart,kend,nkgrid) ;
       
       for thisk=1:nkgrid
           Z(:,thisk) = rfuncirr(avec,kgrid(thisk),theta,P,kstart,kend);
       end
       %contour(meshgrid(kgrid_log*100),meshgrid(theta*100),Z)
       plot((kgrid/kss-1)*100,log10(abs(Z(26,:)/css)))
    end
    
    
    
    %%
    
    
    nreps = 500;
    welfrep=zeros(nreps,1);
    nperiods = 400;
    
    %kgrid_log = -0.05:0.01:0.05;
    
    %kgrid_log = [0.026, 0, 0.10] ;
    %a_start_pos_list = [9, 26, 37];
    
    kgrid_log = [0] ;
    a_start_pos_list = [26];
    
    
    kgrid = exp(kgrid_log+log(kss));
    
    nkgrid = length(kgrid);
    
    
    if compute_welfare
        
        for kindx = 1:nkgrid;
            
            
            for thisrep = 1:nreps
                [kindx thisrep]
                
                %if mod(thisrep,2)==1
                    s=RandStream('mt19937ar','seed',thisrep);
                    RandStream.setDefaultStream(s);
                    %chain=markov_luca(P,nperiods,ceil(nstates/2));
                    chain=markov_luca(P,nperiods,a_start_pos_list(kindx));
           
                %else
                %    chain = nstates+1-chain;
                %end
                k_path = zeros(nperiods+1,1);
                c_path = zeros(nperiods+1,1);
                lambdak_path = zeros(nperiods+1,1);
                resid_path = zeros(nperiods+1,1);
                
                k_path(1) = kgrid(kindx);
                lambdak_path(1) = 0;
                for period = 2:nperiods+1;
                    
                    a_state = [A_sol(:,chain(period-1))]';
                    
                    efunc= a_state*chebypol(k_path(period-1),order,kstart,kend);
                    c_path(period-1) = efunc^(-1/GAMMAC);
                    k_path(period) = invcfunc(theta(chain(period-1)),k_path(period-1),c_path(period-1));
                    
                    if (k_path(period)-(1-DELTAK)*k_path(period-1))<PHI*iss
                        k_path(period) = PHI*iss+(1-DELTAK)*k_path(period-1);
                        c_path(period-1) = cfunc(theta(chain(period-1)),k_path(period-1),k_path(period));
                        
                        
                    end
                    
                end
                
                welfrep(thisrep) = (BETA.^(0:nperiods-1))*(c_path(1:end-1).^(1-GAMMAC))/(1-GAMMAC);
            end
            
            welf(kindx) = mean(welfrep);
            
        end
        save(file_name_welf,'welf')
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Plot decision rules
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    
    if compute_decision_rules
        figure
        kgrid = linspace(kss*.97,kss*(1.3),1000);
        npoints = length(kgrid);
        
        for ia=[ 1:nstates]
            a_state = [A_sol(:,ia)]';
            
            for i  = 1:npoints
                efunc = a_state*chebypol(kgrid(i),order,kstart,kend);
                
                
                kprime(i) = invcfunc(theta(ia),kgrid(i),efunc^(-1/GAMMAC));
                % next verify guess for lambdak
                if kprime(i) -(1-DELTAK)*kgrid(i)<PHI*iss
                    kprime(i) = (1-DELTAK)*kgrid(i)+PHI*iss;
                    
                end
            end
            
            %h=plot((kgrid-kss)/kss*100,(kprime-kgrid)/kss*100); hold on;
            h=plot((kgrid-kss)/kss*100,(kprime-(1-DELTAK)*kgrid)/iss); hold on;
            
            set(h,'Color',[ia/nstates 0.5 1 ])
        end
        ylabel('Percent Dev. from Previous Period')
        xlabel('Percent Dev. from SS')
        title('Decision Rules')
    end

end
