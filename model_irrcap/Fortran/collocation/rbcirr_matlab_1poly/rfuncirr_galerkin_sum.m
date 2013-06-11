
function integral = rfuncirr_galerkin_sum(avec,galerkin_nodes,theta,P,kstart,kend)

integral = rfuncirr_galerkin(avec,galerkin_nodes,theta,P,kstart,kend);

integral = 1000*sum(integral.^2);
