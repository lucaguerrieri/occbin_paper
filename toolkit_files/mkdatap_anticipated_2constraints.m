function [zdata]=mkdatap_anticipated_2constraints(nperiods,decrulea,decruleb,...
    cof,Jbarmat,...
    cof10,Jbarmat10,Dbarmat10,...
    cof01,Jbarmat01,Dbarmat01,...
    cof11,Jbarmat11,Dbarmat11,...
    regime1,regimestart1,...
    regime2,regimestart2,...
    violvecbool,endog_,exog_,...
    irfshock,scalefactormod,init)


nvars = size(endog_,1);


if nargin<16
    init=zeros(nvars,1);
end

if nargin<15;
    scalefactormod=1;
end


nshocks = size(irfshock,1);
for i = 1:nshocks
    shockpos = strmatch(irfshock(i,:),exog_,'exact');
    if ~isempty(shockpos)
        irfshockpos(i) = shockpos;
    else
        error(['Shock ',irfshock(i,:),' is not in the model']);
    end
end



Cbarmat = cof(:,1:nvars);
Bbarmat = cof(:,nvars+1:2*nvars);
Abarmat = cof(:,2*nvars+1:3*nvars);


% cofstar contains the system for the model when the constraint binds


Cbarmat10 = cof10(:,1:nvars);
Bbarmat10 = cof10(:,nvars+1:2*nvars);
Abarmat10 = cof10(:,2*nvars+1:3*nvars);

Cbarmat01 = cof01(:,1:nvars);
Bbarmat01 = cof01(:,nvars+1:2*nvars);
Abarmat01 = cof01(:,2*nvars+1:3*nvars);

Cbarmat11 = cof11(:,1:nvars);
Bbarmat11 = cof11(:,nvars+1:2*nvars);
Abarmat11 = cof11(:,2*nvars+1:3*nvars);

% get the time-dependent decision rules
nregimes1 = length(regime1);
nregimes2 = length(regime2);

Tmax = max([regimestart1(nregimes1) regimestart2(nregimes2)])-1;  % Tmax is the position of the last period
% when the constraint binds

if Tmax > 0
    P = zeros(nvars,nvars,Tmax);
    D = zeros(nvars,Tmax);
    
%     invmat = inv((Astarbarmat*decrulea+Bstarbarmat));
%     P(:,:,Tmax) = -invmat*Cstarbarmat;
%     D(:,Tmax) = -invmat*Dstarbarmat;
    
    
    if (violvecbool(Tmax,1) & ~violvecbool(Tmax,2))
    %XXX fix next three lines
    invmat = inv((Abarmat10*decrulea+Bbarmat10));
    P(:,:,Tmax) = -invmat*Cbarmat10;
    D(:,Tmax) = -invmat*Dbarmat10;  
    elseif (violvecbool(Tmax,1) & violvecbool(Tmax,2))
    invmat = inv((Abarmat11*decrulea+Bbarmat11));
    P(:,:,Tmax) = -invmat*Cbarmat11;
    D(:,Tmax) = -invmat*Dbarmat11;
    else
    invmat = inv((Abarmat01*decrulea+Bbarmat01));
    P(:,:,Tmax) = -invmat*Cbarmat01;
    D(:,Tmax) = -invmat*Dbarmat01;  
    end
    
    invmat = inv(Bbarmat);
    Amat = -invmat*Abarmat;
    Cmat = -invmat*Cbarmat;
    Jmat = -invmat*Jbarmat;
    
    invmat = inv(Bbarmat01);
    Amat01 = -invmat*Abarmat01;
    Cmat01 = -invmat*Cbarmat01;
    Dmat01 = -invmat*Dbarmat01;
    Jmat01 = -invmat*Jbarmat01;

    invmat = inv(Bbarmat11);
    Amat11 = -invmat*Abarmat11;
    Cmat11 = -invmat*Cbarmat11;
    Dmat11 = -invmat*Dbarmat11;
    Jmat11 = -invmat*Jbarmat11;    
    
    invmat = inv(Bbarmat10);
    Amat10 = -invmat*Abarmat10;
    Cmat10 = -invmat*Cbarmat10;
    Dmat10 = -invmat*Dbarmat10;
    Jmat10 = -invmat*Jbarmat10;

    
    
    for i = Tmax-1:-1:1        
        
        if (violvecbool(i,1) & ~violvecbool(i,2))
            invmat = inv(eye(nvars)-Amat10*P(:,:,i+1));
            P(:,:,i)=invmat*Cmat10;
            D(:,i) = invmat*(Amat10*D(:,i+1)+Dmat10);
        elseif (~violvecbool(i,1) & violvecbool(i,2))
            invmat = inv(eye(nvars)-Amat01*P(:,:,i+1));
            P(:,:,i)=invmat*Cmat01;
            D(:,i) = invmat*(Amat01*D(:,i+1)+Dmat01);
        elseif (violvecbool(i,1) & violvecbool(i,2))
            invmat = inv(eye(nvars)-Amat11*P(:,:,i+1));
            P(:,:,i)=invmat*Cmat11;
            D(:,i) = invmat*(Amat11*D(:,i+1)+Dmat11);
        else
            invmat = inv(eye(nvars)-Amat*P(:,:,i+1));
            P(:,:,i)=invmat*Cmat;
            D(:,i) = invmat*(Amat*D(:,i+1));
        end
        
    end

    
% Double check the appropriate invmat in each case
% right now -- inherited from previous loop
if Tmax > 1
    
    if ( ~violvecbool(1,1) & violvecbool(1,2) )
        E = invmat*Jmat01;
    elseif ( violvecbool(1,1) & ~violvecbool(1,2) )
        E = invmat*Jmat10;
    elseif ( violvecbool(1,1) & violvecbool(1,2) )
        E = invmat*Jmat11;
    else
        E = invmat*Jmat;
    end
    
else  % Tmax is equal to 1    
%     invmat = inv((Astarbarmat*decrulea+Bstarbarmat));
%     E = -invmat*Jstarbarmat;
    
    if ( ~violvecbool(1,1) & violvecbool(1,2) )
       invmat = inv((Abarmat01*decrulea+Bbarmat01));
       E = -invmat*Jbarmat01;  
    elseif ( violvecbool(1,1) & violvecbool(1,2) )
        invmat = inv((Abarmat11*decrulea+Bbarmat11));
        E = -invmat*Jbarmat11;     
    else
        invmat = inv((Abarmat10*decrulea+Bbarmat10));
        E = -invmat*Jbarmat10;
            
    end
    
end

    
end

% generate data
% history will contain data, the state vector at each period in time will
% be stored columnwise.
history = zeros(nvars,nperiods+1);
history(:,1) = init;
errvec = zeros(size(exog_,1),1);

for i = 1:nshocks
    errvec(irfshockpos(i)) = scalefactormod(i);
end

% deal with shocks
irfpos =1;
if irfpos <=Tmax
    history(:,irfpos+1) = P(:,:,irfpos)* history(:,irfpos)+...
        D(:,irfpos) + E*errvec;
else
    history(:,irfpos+1) = decrulea*history(:,irfpos)+decruleb*errvec;
end

% all other periods
for irfpos=2:nperiods+1
    if irfpos <=Tmax
        history(:,irfpos+1) = P(:,:,irfpos)* history(:,irfpos)+...
            D(:,irfpos);
    else
        history(:,irfpos+1) = decrulea*history(:,irfpos);
    end
end


history=history';
zdata = history(2:end,:);