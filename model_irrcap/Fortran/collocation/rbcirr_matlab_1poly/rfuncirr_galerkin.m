
function integral = rfuncirr_galerkin(avec,galerkin_nodes,theta,P,kstart,kend)

% initialize

nstates = length(theta);
nparams = length(avec)/nstates;
integral = zeros(nparams*nstates,1);
nnodes = length(galerkin_nodes);


Rvec = rfuncirr2(avec,galerkin_nodes,theta,P,kstart,kend);
R = reshape(Rvec,nnodes,nstates);


chebymat = zeros(nparams+1,nnodes);
order = nparams;
for i = 1:nnodes
    chebymat(:,i) = chebypol(galerkin_nodes(i),order,kstart,kend);
end
chebymat = chebymat(2:end,:);

%for i = 1:nparams
%   integral((i-1)*nstates+(1:nstates))=transpose(sum(reshape(repmat(chebymat(i,:)',nstates,1).*Rvec,nnodes,nstates))); 
%end

j = 0; 
for k = 1:nstates
for i = 1:nparams
  
   j = j+1;
   %integral((i-1)*nstates+(1:nstates))=transpose(sum(reshape(repmat(chebymat(i,:)',nstates,1).*Rvec,nnodes,nstates))); 
   integral(j) = R(:,k)'*chebymat(i,:)';
   
end
end

integral=integral*100;