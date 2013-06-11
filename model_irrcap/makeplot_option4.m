  
  % put things back in levels
  % _l denotes linearized solution
  % _p denotes piece-wise linear
  % _s denotes the nonlinear stochastic solution
  % _d denotes the nonlinear perfect foresight solution
  
  disp(' ')
  disp(' ')
  disp(' ')
  disp('Comparison with FULL NON-LINEAR MODEL')
  
  % Create levels
  c_l=exp(c_uncdifference+log(css));
  c_p=exp(c_difference+log(css));
  c_s=simC';
 
  i_l=exp(i_uncdifference+log(iss));
  i_p=exp(i_difference+log(iss));
  i_s=simI';
 
  
  
  S=1;
  fmt='%6.4f';
  
  disp(' ')
  disp('Investment')
  disp('METHOD     mean    st.dev.  skewness   ')
  disp(['Linear      ' num2str(mean((i_l)),fmt) '  ' num2str(S*std(log(i_l)),fmt) '  ' num2str(skewness(log(i_l)),fmt) '  '])
  disp(['PLinear     ' num2str(mean((i_p)),fmt) '  ' num2str(S*std(log(i_p)),fmt) '  ' num2str(skewness(log(i_p)),fmt) '  '])
  %disp(['NLinear Det ' num2str(mean((i_d)),fmt) '  ' num2str(S*std(log(i_d)),fmt) '  ' num2str(skewness(log(i_d)),fmt) '  '])
  disp(['NLinear Sto ' num2str(mean((i_s)),fmt) '  ' num2str(S*std(log(i_s)),fmt) '  ' num2str(skewness(log(i_s)),fmt) '  '])
  
  disp(' ')
  disp('Consumption')
  disp('METHOD     mean    st.dev.  skewness   ')
  disp(['Linear      ' num2str(mean((c_l)),fmt) '  ' num2str(S*std(log(c_l)),fmt) '  ' num2str(skewness(log(c_l)),fmt) '  '])
  disp(['PLinear     ' num2str(mean((c_p)),fmt) '  ' num2str(S*std(log(c_p)),fmt) '  ' num2str(skewness(log(c_p)),fmt) '  '])
  %disp(['NLinear Det ' num2str(mean((c_d)),fmt) '  ' num2str(S*std(log(c_d)),fmt) '  ' num2str(skewness(log(c_d)),fmt) '  '])
  disp(['NLinear Sto ' num2str(mean((c_s)),fmt) '  ' num2str(S*std(log(c_s)),fmt) '  ' num2str(skewness(log(c_s)),fmt) '  '])

 
  

   disp(' ')
   disp('Correlation between investment and consumption ')
  cor_l=corrcoef(log(i_l),log(c_l));
  cor_p=corrcoef(log(i_p),log(c_p));
  cor_s=corrcoef(log(i_s),log(c_s));
  disp(['Linear      ' num2str(cor_l(1,2),fmt)  ])
  disp(['PLinear     ' num2str(cor_p(1,2),fmt)  ])
  disp(['NLinear Sto ' num2str(cor_s(1,2),fmt)  ])

  
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%   % plot cdf 
%   figure(100)
%   cdf_p=cdfplot(c_p); 
%   c_p_xdata = get(cdf_p,'XData');
%   c_p_ydata = get(cdf_p,'YData');
%   cdf_l=cdfplot(c_l); 
%   c_l_xdata = get(cdf_l,'XData');
%   c_l_ydata = get(cdf_l,'YData');
%   cdf_s=cdfplot(c_s); 
%   c_s_xdata = get(cdf_s,'XData');
%   c_s_ydata = get(cdf_s,'YData');
%   cdf_d=cdfplot(c_d); 
%   c_d_xdata = get(cdf_d,'XData');
%   c_d_ydata = get(cdf_d,'YData');
% 
%   
%   close(100)
% 


[ns xs] = get_pdf((i_s/iss-1)*100,30);
[np xp] = get_pdf((i_p/iss-1)*100,30);
[nl xl] = get_pdf((i_l/iss-1)*100,30);


subplot(2,1,1)

plot(xp,np*(xp(2)-xp(1)),'k','LineWidth',2); hold on
plot(xs,ns*(xs(2)-xs(1)),'r--','LineWidth',2); hold on
plot(xl,nl*(xl(2)-xl(1)),'b-.','LineWidth',2); hold off
ylim([0, 1])
title('Investment: Cumulative Distribution Function')
xlabel('Investment, percent deviation from steady state')
legend(legendlist)
grid on


[ns xs] = get_pdf((c_s/css-1)*100,30);
[np xp] = get_pdf((c_p/css-1)*100,30);
[nl xl] = get_pdf((c_l/css-1)*100,30);


subplot(2,1,2)

plot(xp,np*(xp(2)-xp(1)),'k','LineWidth',2); hold on
plot(xs,ns*(xs(2)-xs(1)),'r--','LineWidth',2); hold on
plot(xl,nl*(xl(2)-xl(1)),'b-.','LineWidth',2); hold off

ylim([0, 1])
title('Consumption: Cumulative Distribution Function')
xlabel('Consumtion, percent deviation from steady state')
grid on

%   figure
%   plot(c_p_xdata,c_p_ydata,'k','Linewidth',2); hold on
%   plot(c_d_xdata,c_d_ydata,'b--','Linewidth',2); hold on
%   plot(c_s_xdata,c_s_ydata,'r-.','Linewidth',2); hold on
%   plot(c_l_xdata,c_l_ydata,'m','Linewidth',2); hold on
%   xlabel('Consumption')
%   ylabel('Cumulative Distribution Function of Consumption')
%   legend('Piecewise Linear','Nonlinear Perfect Foresight','Nonlinear Stochastic','Linearized');
% 
%   
%   % run Kolmogorov Smirnov test
%   [h_ps,pval_ps] = kstest2(c_p,c_s)
%   [h_ls,pval_ls] = kstest2(c_l,c_s)
%   [h_pd,pval_pd] = kstest2(c_p,c_d)
%   [h_ld,pval_ld] = kstest2(c_l,c_d)
%   
%   % Now compute average value function using stochastic value function as benchmark
%   addpath('C:\E\occbin_beta\nonlinear_models\houseprice')
%   if nh==1
%   v_l=interpn(Q,B,V,q_l,lg(b_l,1)) ;
%   v_p=interpn(Q,B,V,q_p,lg(b_p,1)) ;
%   v_d=interpn(Q,B,V,q_d,lg(b_d,1)) ;
%   v_s=interpn(Q,B,V,q_s,lg(b_s,1)) ;
%   else
%   v_l=interpn(Q,B,H,V,q_l,lg(b_l,1),lg(h_l,1)) ;
%   v_p=interpn(Q,B,H,V,q_p,lg(b_p,1),lg(h_p,1)) ;
%   v_d=interpn(Q,B,H,V,q_d,lg(b_d,1),lg(h_d,1)) ;
%   v_s=interpn(Q,B,H,V,q_s,lg(b_s,1),lg(h_s,1)) ;    
%   end
%   rmpath('C:\E\occbin_beta\nonlinear_models\houseprice')
%   
%   
%     
% %   % Find extra % change in c needed to match c_s2  
% 
%   options= optimset('TolX',0.000000000000001);
%   f = @(db)dexpectedutility(db,b_l,h_l,q_l, b_s,h_s,q_s, Q,B,H,V,V);
%   db_l=fzero(f,0.001,options)/exp(c_ss)*100*(1-BETA);
%   f = @(db)dexpectedutility(db,b_p,h_p,q_p, b_s,h_s,q_s, Q,B,H,V,V);
%   db_p=fzero(f,0.001,options)/exp(c_ss)*100*(1-BETA);
%   f = @(db)dexpectedutility(db,b_d,h_d,q_d, b_s,h_s,q_s, Q,B,H,V,V);
%   db_d=fzero(f,0.001,options)/exp(c_ss)*100*(1-BETA);
%   f = @(db)dexpectedutility(db,b_s,h_s,q_s, b_s,h_s,q_s, Q,B,H,V,V);
%   db_s=fzero(f,0.001,options)/exp(c_ss)*100*(1-BETA);
% 
%   
%   table_appendix = ...
%     [ nanmean(v_l) db_l std(log(c_l)) skewness(log(c_l)) cor_l cor_l1 mean_lev_l 
%       nanmean(v_p) db_p std(log(c_p)) skewness(log(c_p)) cor_p cor_p1 mean_lev_p 
%       nanmean(v_d) db_d std(log(c_d)) skewness(log(c_d)) cor_d cor_d1 mean_lev_d 
%       nanmean(v_s) db_s std(log(c_s)) skewness(log(c_s)) cor_s cor_s1 mean_lev_s ];
%   
%   
% %   table_appendix = ...
% %     [ mean(v_l) mean((c_l2)) mean(h_l) std(log(c_l)) skewness(log(c_l)) cor_l cor_l1 mean_lev_l 
% %       mean(v_p) mean((c_p2)) mean(h_p) std(log(c_p)) skewness(log(c_p)) cor_p cor_p1 mean_lev_p 
% %       mean(v_d) mean((c_d2)) mean(h_d) std(log(c_d)) skewness(log(c_d)) cor_d cor_d1 mean_lev_d 
% %       mean(v_s) mean((c_s2)) mean(h_s) std(log(c_s)) skewness(log(c_s)) cor_s cor_s1 mean_lev_s ];
%   