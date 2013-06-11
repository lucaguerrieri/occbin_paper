function tsflag = tschk(x)

% Syntax:  tsflag = tschk(x)
% Check to see if object x is a timeseries

tsflag = 0;
if(isstruct(x))
	if(isfield(x,'id'))
		if(strcmp(x.id,'timeseries'))
			tsflag = 1;
		end
	end
end

return
