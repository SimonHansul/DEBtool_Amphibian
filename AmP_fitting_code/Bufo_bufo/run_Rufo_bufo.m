close all; 
global pets 
%addpath(pathdef2());
pets = {'Bufo_bufo'}; 
check_my_pet(pets); 

check_my_pet(pets); 
 
estim_options('default');
estim_options('max_step_number',500);
estim_options('max_fun_evals',5e3);


%% Estimation procedure
estim_options('pars_init_method', 2);
estim_options('results_output',2);
estim_options('method', 'nm');

estim_pars;

estim_options('pars_init_method', 1);
close all;
estim_pars; close all;
estim_pars; close all;
estim_pars; close all;
estim_pars; close all;
estim_pars; close all;
estim_pars; close all;
% estim_pars; close all;
% estim_pars; close all;
% estim_pars; close all;
% estim_pars; close all;
% estim_pars; close all;
% estim_pars; close all;
estim_options('max_step_number',5e3);
estim_options('results_output', 4);
estim_pars; close all;

%% Displaying results
estim_options('method', 'no');
estim_options('pars_init_method', 1);
estim_options('results_output', -4);
estim_pars;
