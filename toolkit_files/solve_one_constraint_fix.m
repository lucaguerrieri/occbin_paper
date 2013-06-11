function [zdata zdataconcatenated oobase_ Mbase_] = ...
         solve_one_constraint(modnam,modnamstar,...
                              constraint, constraint_relax,...
                              shockssequence,irfshock,nperiods,tol0,maxiter)
                          
global M_ oo_                           
                           
errlist = [];                

% solve model
eval(['dynare ',modnam,' noclearall'])
oobase_ = oo_;
Mbase_ = M_;
setss

eval(['dynare ',modnamstar,' noclearall'])
oostar_ = oo_;
Mstar_ = M_;

% check inputs
if ~strcmp(Mbase_.endo_names,Mstar_.endo_names)
    error('The two .mod files need to have exactly the same endogenous variables declared in the same order')
end

if ~strcmp(Mbase_.exo_names,Mstar_.exo_names)
    error('The two .mod files need to have exactly the same exogenous variables declared in the same order')
end

if ~strcmp(Mbase_.param_names,Mstar_.param_names)
    warning('The parameter list does not match across .mod files')
else
    if max(abs(Mbase_.params-Mstar_.params))>tol0
       warning('The two .mod files have different parameter values') 
    end
end



nvars = Mbase_.endo_nbr;
ys_ = oobase_.dr.ys;


[hm1,h,hl1,Jbarmat] = get_deriv(Mbase_,ys_);
cof = [hm1,h,hl1];

[hm1,h,hl1,Jstarbarmat,resid] = get_deriv(Mstar_,ys_);
cofstar = [hm1,h,hl1];
Dstarbarmat = resid;


[decrulea,decruleb]=get_pq(oobase_.dr);
endog_ = M_.endo_names;
exog_ =  M_.exo_names;


% processes the constrain so as to uppend a suffix to each
% endogenous variables
constraint_difference = process_constraint(constraint,'_difference',Mbase_.endo_names,0);

constraint_relax_difference = process_constraint(constraint_relax,'_difference',Mbase_.endo_names,0);




nshocks = size(shockssequence,1);
init = zeros(nvars,1);
zdataconcatenated = zeros(nperiods,nvars);
wishlist = endog_;
  
violvecbool = zeros(nperiods+1,1);
for ishock = 1:nshocks
  % generate IRFs for wishlist. IRFs are save in ZDATA
  zdata = mkdata(nperiods+1,decrulea,decruleb,endog_,exog_,wishlist,irfshock,shockssequence(ishock),init);
  
%   [regime regimestart]=map_regime(violvecbool);
%   
%   [zdata]=mkdatap_anticipated(nperiods,decrulea,decruleb,...
%         cof,Jbarmat,cofstar,Jstarbarmat,Dstarbarmat,...
%         regime,regimestart,violvecbool,...
%         endog_,exog_,irfshock,shockssequence(ishock,:),init);
      
  
  
  % unpack ZDATA
  nwishes = size(wishlist,1);
  for i=1:nwishes
    eval([deblank(wishlist(i,:)),'_difference=zdata(:,i);']);
  end
  
  % find periods in which the constraint is binding
  violvecbool = eval(constraint_difference);
  violvec = find(violvecbool);
  
  
  changes=1;
  iter = 0;
 
  if isempty(violvec)
    changes = 0;
  else
   
    while (changes & iter<maxiter)
      iter = iter +1;
      
      
      [regime regimestart]=map_regime(violvecbool);
      
      
      [zdata]=mkdatap_anticipated(nperiods,decrulea,decruleb,...
        cof,Jbarmat,cofstar,Jstarbarmat,Dstarbarmat,...
        regime,regimestart,violvecbool,...
        endog_,exog_,irfshock,shockssequence(ishock,:),init);
      
      for i=1:nwishes
        eval([deblank(wishlist(i,:)),'_difference=zdata(:,i);']);
      end
      
   
      
      newviolvecbool = eval(constraint_difference);  
      relaxconstraint = eval(constraint_relax_difference);
            
      
      
      % check if changes
      if (max(newviolvecbool-violvecbool>0)) | sum(relaxconstraint(find(violvecbool==1))>0)
        changes = 1;
      else
        changes = 0;
      end
      
      
      violvecbool = (violvecbool|newviolvecbool)-(relaxconstraint & violvecbool);
      
     
    end
  end
  init = zdata(1,:);
  zdataconcatenated(ishock,:)=init;
  init= init';
  
  % reset violvecbool for next period -- consistent with expecting no
  % additional shocks
%   violvecbool=[violvecbool(2:end);0];
  
end

if changes ==1
  display('Did not converge -- increase maxiter')
end

zdataconcatenated(ishock+1:end,:)=zdata(2:nperiods-ishock+1,:);

zdata = mkdata(nperiods,decrulea,decruleb,endog_,exog_,wishlist,irfshock,shockssequence);

