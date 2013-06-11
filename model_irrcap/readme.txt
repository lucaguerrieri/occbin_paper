The original programs are in /Users/Jason/Documents/MATLAB/consumption/occbin_work/occbin_work

This directory contains a subset of the final programs for replication purposes.




Figures

NB:  For Figures 1 to 6 cd into model_irrcap
     For Figures 7 and 8 cd into XXX

Figure 1   Matteo

Figure 2   IRF to up and down shock
           cd fortran/dp/matlab

           setup irrcap.m with the following options

           model_version = 4;    
           compute_welfare = 0;  
           makefigures = 1;         
           mode = 1;        

           run irrcap.m (if a new fortran solution is desired, modify 
           fortran/dp/sargent_fortran/runsim_rbc.f90, but remember to line
           up the options and parameters across programs, and to update the path where
           the decision rule is saved by the fortran program).

           run runsim_irrcap_new to superimpose the PL solution 
           over the nonlinear solution.

Figure 3   CDF
           cd fortran/dp/matlab

           setup irrcap.m with the following options

           model_version = 4;    
           compute_welfare = 0;  
           makefigures = 1;         
           mode = 2;        

           run irrcap.m (if a new fortran solution is desired, modify 
           fortran/dp/sargent_fortran/runsim_rbc.f90, but remember to line
           up the options and parameters across programs, and to update the path where
           the decision rule is saved by the fortran program).

           Back in the model_irrcap directory, 
           run runsim_irrcap_new.m to superimpose the PL solution 
           over the nonlinear solution.

Figure 4   plot_euler_errors_pl_tech_levels.m
           

Figure 5   cd Fortran/collocation/rbcirr_matlab_1poly
           plot_euler_errors.m

Figure 6   cd Fortran/collocation/rbcirr_matlab_1poly

           run runsim_rbciirr2.m  with the following settings:
           generate_data = 1;
           option = 1;  
           compute_euler_residuals = 1;
  
           compute_welfare=0;
           euler_error_mesh = 0;
           compute_decision_rules = 0;
           load_from_fortran_program = 1;
           version = 1;

           The program loads the decision rule saved by 
           Fortran/collocation/rbcirr_mm/runsim_rbcirr.f90
           (remember to sync output directory and options across fortran
            and matlab, in case a new decision rule is desired).

Figure 7   Matteo
    
Figure 8   Matteo  



Tables:

NB: for tables 1 to 3 cd into model_irrcap
    for tables 6 and 7 cd into XXX


Table 1    Calibration

Table 2    Moments -- cd fortran/dp/matlab

           run irrcap.m with the following settings:

           model_version = 4;    
           compute_welfare = 0;  
           makefigures = 1;         
           mode = 2;        

Table 3   
  
Table 4    Calibration

Table 5

Table 6




           
