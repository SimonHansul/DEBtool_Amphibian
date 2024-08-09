function [par, metaPar, txtPar] = pars_init_Rana_temporaria(metaData)

metaPar.model = 'std'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.z = 1.431;        free.z     = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.01691;      free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.6293;     free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 1052;       free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 7310;       free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hh = 0.09867;   free.E_Hh  = 1;   units.E_Hh = 'J';         label.E_Hh = 'maturity at hatch'; 
par.E_Hb = 7.095;     free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_H42 = 37.09;      free.E_H42  = 1;   units.E_H42 = 'J';         label.E_H42 = 'maturity at metam begin';
par.E_Hj = 79.74;       free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam end';
par.E_Hp = 5.255e+05; free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 3.986e-09;  free.h_a   = 1;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
 
par.T_A = 5455;       free.T_A   = 1;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.TrefL = 273.15;      free.TrefL = 0; units.TrefL = 'K';   label.TrefL='ref. temp. Arrhenius for low temperatures';
par.TA_l = 30000;    free.TA_l = 0;    units.TA_l = 'K';   label.TA_l='Arrhenius temperature low temp.';
par.TrefH = 350.8;      free.TrefH = 0; units.TrefH = 'K';   label.TrefH='ref. temp. Arrhenius for high temperatures';
par.TA_h = 2.e+04 	;    free.TA_h = 0;    units.TA_h = 'K';   label.TA_h='Arrhenius temperature high temp.';
par.del_Mt = 0.04438;     free.del_Mt = 1;   units.del_Mt = '-';        label.del_Mt = 'shape coefficient tadpole'; 
par.del_M = 0.128;      free.del_M = 1;   units.del_M = '-';        label.del_M = 'shape coefficient for frog';
par.del_Msvl=0.1196;     free.del_Msvl = 1;   units.del_Msvl = '-';  label.del_Msvl = 'shape coefficient SVL';
par.f = 1;            free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
par.fAv = 1.292;          free.fAv = 1;   units.fAv = '-';        label.fAv = 'scaled functional response for Aviles data';
par.f_LdL = 0.7303;     free.f_LdL = 1;   units.f_LdL = '-';        label.f_LdL = 'scaled functional response for LdL data'; 
par.fRu = 1.;         free.fRu = 0;   units.fRu = '-';        label.fRu = 'scaled functional response for Ru data';

% par.fL1 = 0.7303;        free.fL1  = 1;   units.fL1 = '-';         label.fL1 = 'scaled functional response for Laurila data low food'; 
% par.fL2 = 1.;        free.fL2  = 0;   units.fL2 = '-';         label.fL2 = 'scaled functional response for Laurila data high food'; 

par.bn = 40.68; free.bn=1;units.bn='-';label.bn='scaling factor for natural temp. stress';
% par.shiftRU = 1;free.shiftRU = 1; units.shiftRU='-'; label.shiftRU='shift time axis';
% par.shiftAV = 1;free.shiftAV = 1; units.shiftAV='-'; label.shiftRU='shift time axis';
%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 


%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
