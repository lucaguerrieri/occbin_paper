function [state_lo state_hi state_lo_weight state_hi_weight] = find_weights(this_theta, theta);


state_lo = max(find(theta<this_theta));
state_hi = min(find(theta>this_theta));

if isempty(state_lo)
    state_lo = state_hi;
    state_lo_weight = 0;
    state_hi_weight = 1;
elseif isempty(state_hi)
    state_hi = state_lo;
    state_hi_weight = 0;
    state_lo_weight = 1;
else
    state_hi_weight = (this_theta-theta(state_lo))/(theta(state_hi)-theta(state_lo));
    state_lo_weight = (theta(state_hi)-this_theta)/(theta(state_hi)-theta(state_lo)); 
end