% quarter_to_annual(q,1): Transform a series q from quarter to annual, using averages
% quarter_to_annual(q,2): Transform a series q from quarter to annual, using end of period
% If series has # entries not multiple of 4, it is assumed that first obs is 1st quarter
% Does no input checking

function a = quarter_to_annual(q,x)

% Monthly data
% disp('Size of q is');
% disp(numel(q));

for iq=4:numel(q)
    if x==1;    aq(iq-3) = (q(iq)+q(iq-1)+q(iq-2)+q(iq-3))/4; end
    if x==2;    aq(iq-3) = q(iq); end
end


iq=1;
for j=1:4:numel(aq)
    a(iq) = aq(j);
    iq=iq+1;
end

% Transform a so that it is a column
a=a';

% disp('Size of a is')
% disp(numel(a))
% 
% plot(a)