
function Rvec = rfuncirr_sum(avec,nodes1,nodes2,theta,P,kstart,kend)

Rvec = rfuncirr(avec,nodes1,nodes2,theta,P,kstart,kend);

Rvec = sum(Rvec.^2);



