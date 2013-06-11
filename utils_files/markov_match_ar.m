function [chain]=markov_match_ar(rho,epsilon,theta,nperiods);


path_ar = zeros(nperiods,1);

path_ar(1) = epsilon;
for i=2:nperiods
   path_ar(i)=rho*path_ar(i-1); 
end


for i=1:nperiods
    if path_ar(i) >= 0;
        chain(i)=max(find(min(abs(path_ar(i)-theta))==abs(path_ar(i)-theta) ));
    else
        chain(i)=min(find(min(abs(path_ar(i)-theta))==abs(path_ar(i)-theta) ));
    end
    
end

