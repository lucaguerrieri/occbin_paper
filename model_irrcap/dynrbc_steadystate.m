function [ys,check] = dynrbc_steadystate(junk,ys)


global M_  


check = 0;

paramfile


def_parm

% export parameters
nparams = size(M_.param_names,1);
for icount = 1:nparams
    eval(['M_.params(icount) = ',M_.param_names(icount,:),';'])
end

% transfer parameters to Dynare.
% send chip to the parameter list
%M_.params(strmatch('chip',M_.param_names,'exact')) = chip;
%M_.params(strmatch('pbss',M_.param_names,'exact')) = pbss;
%M_.params(strmatch('puss',M_.param_names,'exact')) = puss;
%M_.params(strmatch('pcss',M_.param_names,'exact')) = pcss;

klag = log(kss)
k = log(kss);
c = log(css);
i = log(iss);
lambdak = 0;
elambdak = 0;
a=0;
v = log(vss);



nvars = M_.endo_nbr;

ys = zeros(nvars,1);

for i_indx = 1:nvars
   eval(['ys(i_indx)=',M_.endo_names(i_indx,:),';']) 
end


