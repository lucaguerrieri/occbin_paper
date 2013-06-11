%% set inputs for solution below
%  The program produces responses to the shocks selected below under
%  irfshock. Paths for all the endogenous variables in the model selected
%  are produced. VariableName_difference holds the piece-wise linear
%  solution for
%  VariableName.  VariableName_uncdifference holds the linear solution for
%  VariableName.

clear
setpathdynare4

global M_ oo_
global BETA ALPHA RHOA DELTAK GAMMAC PSI PHI kss iss

% modname below chooses model
% directory. But simple param choices are made from paramfile in current
% directory.
modnam = 'dynrbc';
modnamstar = 'dynrbcirr';


compute_euler_residuals = 0;

% needed to save linear decision rule (used to initialize DP solution)
%fortran_path = '/Users/Jason/Documents/MATLAB/consumption/occbin_20121229/model_irrcap/fortran/dp/sargent_fortran/';
fortran_path = 'G:\ofs\research2\Guerrieri\sargent_fortran3\';
save_initial_condition = 0;   % set to one to save the linear decision rule for the model
% to be used to initialize the nonlinear rule
% for the dynamic programming solution


if ispc
    !copy paramfile_gammac_2.m paramfile.m
else
    !cp paramfile_gammac_2.m paramfile.m
end
paramfile
def_parm
SIGMA_EPS = 0.013;
% express the occasionally binding constraint
% in linearized form
% one can use any combination of endogenous variables and parameters
% declared in the the dynare .mod files
% constraint1 defines the first constraint
% if constraint1 is true, solution switches to model2
% but if constraint_relax is true, solution reverts to model1


constraint = 'i<log(PHI)';
constraint_relax ='lambdak<0';

% Pick innovation for IRFs
irfshock =char('erra');      % label for innovation for IRFs
% needs to be an exogenous variable in the
% dynare .mod files

maxiter = 10;
tol0 = 1e-8;


% Option=1: impulse responses
% Option=2: random simulation
% Option=3: comparison with nonlinear model when available

option=4;

%%%%%%%%%%%%%%%% Inputs stop here %%%%%%%%%%%%%%%%%%%%%
%%

if option==1
    nper=1;
    
    shockssequence = [
        zeros(9,1)
        -0.02*ones(nper,1)/nper
        zeros(39,1)
        0.02*ones(nper,1)/nper
        ];         % scale factor for simulations
    nperiods = size(shockssequence,1);            %length of IRFs
    
    
end

if option==2
    nperiods = 25;
    randn('seed',3);
    shockssequence = 1*randn(nperiods,1)*0.03 ;
end


if option==3
    
    % NON-LINEAR MODEL FROM MATTEO
    % =0 Arbitrary sequence T=2000
    % =1 IRF like with extreme values T=100
    % =2 No shocks, just check
    % =3 IRF like with less extreme values
    
    load C:\E\Consumption\borrcon_nonlin\irrcap2;
    overwrite=0; T=100;
    run('C:\E\Consumption\borrcon_nonlin\irrcap1_sim');
    
    save irrcap1_play
    
    EPSILON(1)=log(simZ(1));
    for t=2:numel(simZ)
        EPSILON(t)=(simZ(t)-1)-RHOA*(simZ(t-1)-1);
        EPSILON(t)=log(simZ(t))-RHOA*log(simZ(t-1));
    end
    
    shockssequence = EPSILON';
    nperiods = size(shockssequence,1);
    
    
end

if option==4
    
    % NON-LINEAR MODEL FROM LUCA
    
    load luca_nl_sim;
    
    EPSILON(1)=log(simZ(1));
    for t=2:numel(simZ)
        RHOA = 0.9;
        %EPSILON(t)=(simZ(t)-1)-RHOA*(simZ(t-1)-1);
        % hardwired -- double check for consistency
        EPSILON(t)=log(simZ(t))-RHOA*log(simZ(t-1));
    end
    
    shockssequence = EPSILON';
    nperiods = size(shockssequence,1);
    
    
end


% Solve model, generate model IRFs
[zdata zdataconcatenated oobase_ Mbase_ oostar_ Mstar_] = ...
    solve_one_constraint_temp1(modnam,modnamstar,...
    constraint, constraint_relax,...
    shockssequence,irfshock,nperiods,tol0,maxiter);




if compute_euler_residuals
    kap_difference = zeros(nperiods,1);
    Rvec = zeros(nperiods,1);
    nstates = 51;
    [P, theta] =markovappr(RHOA,SIGMA_EPS,3,nstates);
    apos = strmatch('a',M_.endo_names);
    kpos = strmatch('k',M_.endo_names);
    
    init = 0*zdata(1,:)';
    for thisperiod = 1:nperiods
       
        [Rvec(thisperiod) kap_difference(thisperiod) ] = ...
            rfuncirr_i(init,shockssequence(thisperiod),irfshock,theta,...
                       P,modnam,modnamstar,constraint,constraint_relax,...
                       oobase_,Mbase_,oostar_,Mstar_);
        init = zdataconcatenated(thisperiod,:)';
        init(kpos) = kap_difference(thisperiod);
    end
    
    
    %%% repeat for unconstrained model
    % recycle code above, but setting PHI = 0 -- NB value in Mbase is
    % copied in Mstar by solve1constraint.
    Mbase_.params(strmatch('PHI',Mbase_.param_names,'exact'))=0;
    Rvec_unconstrained = zeros(nperiods,1);
    
    init = 0*zdata(1,:)';
    for thisperiod = 1:nperiods
       
        Rvec_unconstrained(thisperiod) = ...
            rfuncirr_i(init,shockssequence(thisperiod),irfshock,theta,...
                       P,modnam,modnamstar,constraint,constraint_relax,...
                       oobase_,Mbase_,oostar_,Mstar_);
        init = zdata(thisperiod,:)';
        
    end
    
    
end


if save_initial_condition
    
    % save output
    [decrulea,decruleb]=get_pq(oobase_.dr);
    
    endog_ = M_.endo_names;
    exog_ =  M_.exo_names;
    kpos = strmatch('k',endog_)
    apos = strmatch('a',endog_)
    vpos = strmatch('v',endog_)
    nlinvars = size(endog_,1)
    exportmat(decrulea,[fortran_path,'decrulea.dat'])
    exportmat(decruleb,[fortran_path,'decruleb.dat'])
    
    
    eval(['save ',fortran_path,'linsol.mat decrulea decruleb endog_ exog_'])
end

% unpack the IRFs
for i=1:M_.endo_nbr
    eval([deblank(M_.endo_names(i,:)),'_uncdifference=zdata(:,i);']);
    eval([deblank(M_.endo_names(i,:)),'_difference=zdataconcatenated(:,i);']);
end


%% Modify to plot IRFs and decision rules

% modify to plot IRFs
lbss=0;

titlelist = char('Technology','Investment','Consumption','Capital');
percent = 'Percent dev. from steady state';
value = 'value';
ylabels = char(percent,percent,percent,percent);

%i_difference = log(exp(k_difference+log(kss))-(1-DELTAK)*exp(klag_difference+log(kss)))-log(iss);
%i_uncdifference = log(exp(k_uncdifference+log(kss))-(1-DELTAK)*exp(klag_uncdifference+log(kss)))-log(iss);


figtitle = '';
line1=100*[a_difference,i_difference,...
    c_difference,k_difference];
line3=100*[a_uncdifference,i_uncdifference,...
    c_uncdifference,k_uncdifference];



% Figure 1
if (option==3)
    legendlist = cellstr(char('Piecewise Linear','Full nonlinear','Linear'));
    load irrcap1_play;
    line2=100*[log(simZ'),log(simI')-log(iss),...
        log(simC')-log(css),log(simK')-log(kss)];
elseif (option == 4)
    legendlist = cellstr(char('Piecewise Linear','Full nonlinear','Linear'));
    if exist('simLAMBDAK')
        line2=100*[log(simZ'),log(simI')-log(iss),...
            log(simC')-log(css),log(simK')-log(kss)];
    else
        line2=100*[log(simZ'), log(simI')-log(iss),...
            log(simC')-log(css),log(simK') - log(kss)];
    end
else
    legendlist = cellstr(char('Piece-Wise Linear','Linear'));
    line2=NaN*line1;
end

makechart9(titlelist,legendlist,figtitle,-1000,ylabels,line1,line2,line3);




if (option == 4 & compute_euler_residuals == 1)
      
    
    titlelist = char('Technology','Investment','Consumption','Lagrange Multiplier','Capital','Magnitude of Euler Residual');
    percent = 'Percent dev. from steady state';
    value = 'Level';
    percent_c = 'Log 10 Scale';
    ylabels = char(percent,percent,percent,value,percent,percent_c);
    
    %i_difference = log(exp(k_difference+log(kss))-(1-DELTAK)*exp(klag_difference+log(kss)))-log(iss);
    %i_uncdifference = log(exp(k_uncdifference+log(kss))-(1-DELTAK)*exp(klag_uncdifference+log(kss)))-log(iss);
    
    
    figtitle = '';
    line1=100*[a_difference,i_difference,...
        c_difference,lambdak_difference/100,...
        k_difference,log10(abs(Rvec/css))/100];
    
    line3 =100*[a_uncdifference,i_uncdifference,...
        c_uncdifference,lambdak_uncdifference/100,...
        k_uncdifference,log10(abs(Rvec_unconstrained/css))/100];
    
    
    line2=100*[log(simZ'),log(simI')-log(iss),...
        log(simC')-log(css),simLAMBDAK'/100,...
        log(simK')-log(kss),log10(abs(resid_path'))/100];
    
    legendlist = cellstr(char('Piece-Wise Linear','Non-Linear','Log-linear of Unconstrained Model'));
    makechart9(titlelist,legendlist,figtitle,-1000,ylabels,line1,line2,line3);
    
end



if (option==3)
    
    c_l=exp(c_uncdifference)*css;
    c_pl=exp(c_difference)*css;
    c_nl=simC;
    k_l=exp(k_uncdifference)*kss;
    k_nl=simK;
    k_pl=exp(k_difference)*kss;
    a_l=exp(a_uncdifference);
    a_pl=exp(a_difference);
    a_nl=simZ;
    i_l=exp(i_uncdifference)*iss;
    i_pl=exp(i_difference)*iss;
    i_nl=simK-(1-DELTAK)*simKlag;
    
    warning off
    disp(' ')
    disp(' ')
    disp(' ')
    disp('Comparison with FULL NON-LINEAR MODEL')
    
    
    disp(' ')
    disp('Log Consumption in various models')
    disp('            mean    st.dev.   kurtosis   ')
    disp(['Linear     ' num2str(mean(log(c_l(2:end))),4)  '  ' num2str(std(log(c_l(2:end))),4) '  ' num2str(skewness(log(c_l(2:end)),4)) '  '])
    disp(['PLinear    ' num2str(mean(log(c_pl(2:end))),4) '  ' num2str(std(log(c_pl(2:end))),4) '  ' num2str(skewness(log(c_pl(2:end)),4)) '  '])
    disp(['NLinear    ' num2str(mean(log(c_nl(2:end))),4) '  ' num2str(std(log(c_nl(2:end))),4) '  ' num2str(skewness(log(c_nl(2:end)),4)) '  '])
    
    
    disp(' ')
    disp('Log Capital in various models')
    disp('            mean    st.dev.   kurtosis   ')
    disp(['Linear     ' num2str(mean(log(k_l(2:end))),5)  '  ' num2str(std(log(k_l(2:end))),4) '  ' num2str(skewness(log(k_l(2:end)),4)) '  '])
    disp(['PLinear    ' num2str(mean(log(k_pl(2:end))),5) '  ' num2str(std(log(k_pl(2:end))),4) '  ' num2str(skewness(log(k_pl(2:end)),4)) '  '])
    disp(['NLinear    ' num2str(mean(log(k_nl(2:end))),5) '  ' num2str(std(log(k_nl(2:end))),4) '  ' num2str(skewness(log(k_nl(2:end)),4)) '  '])
    
    
    disp(' ')
    disp('TFP in various models')
    disp('            mean    st.dev.   kurtosis   ')
    disp(['Linear     ' num2str(mean(log(a_l(2:end))),4)  '  ' num2str(std(log(a_l(2:end))),4) '  ' num2str(skewness(log(a_l(2:end)),4)) '  '])
    disp(['PLinear    ' num2str(mean(log(a_pl(2:end))),4) '  ' num2str(std(log(a_pl(2:end))),4) '  ' num2str(skewness(log(a_pl(2:end)),4)) '  '])
    disp(['NLinear    ' num2str(mean(log(a_nl(2:end))),4) '  ' num2str(std(log(a_nl(2:end))),4) '  ' num2str(skewness(log(a_nl(2:end)),4)) '  '])
    
    disp(' ')
    disp('% of times when constraint is violated (Linear) or binding    ')
    f_l=numel(find(i_l<PHI*iss*1.01))/T*100;
    f_pl=numel(find(i_pl<PHI*iss*1.01))/T*100;
    f_nl=numel(find(i_nl<PHI*iss*1.01))/T*100;
    disp(['Linear     ' num2str(f_l,4)  ])
    disp(['PLinear    ' num2str(f_pl,4)  ])
    disp(['NLinear    ' num2str(f_nl,4)  ])
    
end


if option==4
    figure
    makeplot_option3
end
