function [par, metaPar, txtPar] = pars_init_Dryophytes_versicolor(metaData)

metaPar.model = 'std'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.z = 1.041;        free.z     = 1;   units.z = '-';            label.z = 'zoom factor for toad'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.02247;      free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.9143;     free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 875.1;       free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 7443;       free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hh = 0.01049;    free.E_Hh  = 1;   units.E_Hh = 'J';         label.E_Hh = 'maturity at hatch'; 
par.E_Hb = 0.1884;    free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_H42 = 39.9;      free.E_H42  = 1;   units.E_H42 = 'J';         label.E_H42 = 'maturity at G42'; 
par.E_Hj = 	54.01;      free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at end of met. climax'; 
par.E_Hp = 1.477e+04	;      free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 1.609e-08;  free.h_a   = 1;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
par.T_A = 6981;       free.T_A   = 1;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.del_M = 0.1797;   free.del_M = 1;   units.del_M = '-';        label.del_M = 'shape coefficient for frog'; 
par.f = 1.0;            free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
% par.fH = 1.0;            free.fH   = 1;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
par.fL = 0.6271;       free.fL     = 1;   units.fL = '-';            label.fL = 'scaled functional response for Smith2004 data'; 
par.fS04 = 0.3161;       free.fS04     = 1;   units.fS04 = '-';            label.fS04 = 'scaled functional response for Smith2004 data'; 
par.fB99L = 0.6309;      free.fB99L     = 1;   units.fB99L = '-';            label.fB99L = 'scaled functional response for Beachy1999 data low food'; 
par.fB99H = 1.;     free.fB99H     = 0;   units.fB99H = '-';            label.fB99H = 'scaled functional response for Beachy1999 data high food'; 
par.bn = 0;       free.bn = 0;      units.bn='-'; label.bn='temperature stress scale';

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 


%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
