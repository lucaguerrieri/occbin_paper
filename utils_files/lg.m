function lagv = lg(v,k)

% Syntax:  lagv = mylag(v,k)
%
% "MyLag" a vector or timeseries v by k periods, i.e. shift it right by k periods, 
% padding the beginning with k values corresponding to INITIAL observation.
% Modified version of Fuhrer's LAG, works only for vectors!

% Default for lag argument k is 1

if(nargin==1)
   k = 1;
end

if(k<0)
   error('Positive lag required')
end

if(k==0)
   lagv = v;   
   return
end


[nrows,ncols] = size(v);

% Check to see if v is a timeseries

tsflag = tschk(v);

if(tsflag)

	% For timeseries, lagged series is just same data as
	% original series but starting at the original startdate
	% shifted forward k periods.

	savev = v;
	freq = getfq_ts(savev);
	sd  = getst_ts(savev);
	v = getdt_ts(v);

	% Store back in timeseries matrix with new start date if input is a
	% timeseries matrix 

	sd = tsshift(sd,k,freq);   % New startdate

	kstr = int2str(k);
	lkstr = length(kstr);
	v = tseries(v,sd,freq,['L',kstr,' ',setstr(getnm_ts(savev))]);  

	lagv = v;

else

	% Transpose column vectors for shiftright operation

	if(nrows > ncols)
	  v = v';
    end

	lagv = shftrght(v,k);
    lagv(1:k) = lagv(k+1) ;

	% Put result back into column vector form if original was a column
	% vector 

	if(nrows > ncols)
	  lagv = lagv';
	end

end

