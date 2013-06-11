function integral = rfunc_mm(avec,galerkin_nodes,theta,P,kstart,kend)

% initialize

nstates = length(theta);
nparams = length(avec)/nstates;
integral = zeros(nparams*nstates,1);

Rvec = rfunc(avec,galerkin_nodes,theta,P,kstart,kend);

nnodes = length(galerkin_nodes);
polmat = zeros(nnodes,nparams);
order = nparams;
for i = 1:nparams
    polmat(:,i) = galerkin_nodes.^(i-1);
end


R = reshape(Rvec,nnodes,nstates);

j = 0; 
for k = 1:nstates
for i = 1:nparams
  
   j = j+1;
   %integral((i-1)*nstates+(1:nstates))=transpose(sum(reshape(repmat(chebymat(i,:)',nstates,1).*Rvec,nnodes,nstates))); 
   integral(j) = R(:,k)'*polmat(:,i);
   
end
end

integral=integral;