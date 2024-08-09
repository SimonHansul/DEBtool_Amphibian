function [par, metaPar, txtPar] = pars_init_Xenopus_laevis(metaData)

metaPar.model = 'std'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.T_A = 4948;       free.T_A   = 1;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.TrefL = 273.15;   free.TrefL = 0;   units.TrefL = 'K';        label.TrefL='ref. temp. Arrhenius for low temperatures';
par.TA_l = 30000;     free.TA_l  = 0;   units.TA_l = 'K';         label.TA_l='Arrhenius temperature low temp.';
par.TrefH = 361.15;   free.TrefH = 0;   units.TrefH = 'K';        label.TrefH='ref. temp. Arrhenius for high temperatures';
par.TA_h = 30000;     free.TA_h  = 0;   units.TA_h = 'K';         label.TA_h='Arrhenius temperature high temp.';
par.z = 3.708;         free.z     = 1;   units.z = '-';            label.z = 'zoom factor for females'; 
par.F_m = 6.5;        free.F_m   = 0;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;      free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;      free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.04732;      free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.9391;     free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma';
par.kap_R = 0.95;     free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 209;      free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;          free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;      free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 7385;       free.E_G   = 1;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb =  0.2686;    free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hp = 4.327e+04; free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 1.53e-10;  free.h_a   = 1;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;     free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% other parameters 
% parV = [138.3 181.5 296.7 413.4 497 835.2 1270 1569 1848 1973 2129 2340 2025]; %parameter vector maturity at stages NF50 - NF62
% for i = 1:13  
%    E_Hx = ['E_H',num2str(i+49)];
%    par.(E_Hx) = parV(i);
%    free.(E_Hx)  = 1;   units.(E_Hx) = 'J';    label.(E_Hx) = ['maturity at stage NF',num2str(i+49)];    
% end
%par.E_H62 = 166.9;     free.E_H62  = 1;   units.E_H62 = 'J';         label.E_H62 = 'maturity at start of climax'; 
par.E_H62 = 127.4;     free.E_H62  = 1;   units.E_H62 = 'J';         label.E_H62 = 'maturity at start of climax'; %TEST
par.E_Hh = 0.02497;     free.E_Hh  = 1;   units.E_Hh = 'J';         label.E_Hh = 'maturity at hatch'; 
%par.E_Hj = 241.3;      free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
par.E_Hj = 249;      free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; %TEST
par.del_M = 0.3312;    free.del_M = 1;   units.del_M = '-';        label.del_M = 'shape coefficient for frog'; 
par.del_Me = 0.1753;   free.del_Me = 1;  units.del_Me = '-';       label.del_Me = 'shape coefficient for SVL length for tadpole'; 
par.del_Meh = 0.1503;  free.del_Meh = 1; units.del_Meh = '-';      label.del_Meh = 'shape coefficient for hind limb length of tadpole';
par.del_Mt = 0.05475;  free.del_Mt = 1;  units.del_Mt = '-';       label.del_Mt = 'shape coefficient for tadpole'; 
par.f = 1;             free.f     = 0;   units.f = '-';            label.f = 'scaled functional response for 0-var data'; 
par.fMAR = 1.00;      free.fMAR  = 0;   units.fMAR = '-';            label.fMAR = 'scaled functional response for AMA data of Marini et al. 2023'; 
par.Xkratio = 0.2477;  free.Xkratio = 1; units.Xkratio = '-';      label.Xkratio = 'ratio between food density at half feeding rate and food density';
par.f_tL = 1.03;      free.f_tL  = 1;   units.f_tL = '-';         label.f_tL = 'scaled functional response for tL,tW data'; 
par.f_tad = 0.6951;    free.f_tad = 1;   units.f_tad = '-';        label.f_tad = 'scaled functional response for tadpole data'; 
% par.fWa = 1;            free.fWa = 0; units.fWa = '-'; label.fWa = 'scaled functional response for Walsh08';
par.fRu = 1;           free.fRu = 0; units.fRu = '-'; label.fRu = 'scaled functional response for Ruthsatz22';
par.fNa = 1.0;        free.fNa = 0; units.fNa = '-'; label.fNa = 'scaled functional response for Nations11';
par.shiftNa = 2;       free.shiftNa=0; units.shiftNa = 'd'; label.shiftNa = 'shift of the origin with respect to the hatching for Nations11';
% par.f_tad1 = 0.9453;  free.f_tad1 = 1;  units.f_tad1 = '-';       label.f_tad1 = 'scaled functional response for tadpole data set 1'; 
% par.f_tad2 = 0.9718;  free.f_tad2 = 1;  units.f_tad2 = '-';       label.f_tad2 = 'scaled functional response for tadpole data set 2'; 
% par.f_tad3 = 1.239;   free.f_tad3 = 1;  units.f_tad3 = '-';       label.f_tad3 = 'scaled functional response for tadpole data set 3'; 
% par.f_tad4 = 1.123;   free.f_tad4 = 1;  units.f_tad4 = '-';       label.f_tad4 = 'scaled functional response for tadpole data set 4'; 
% par.f_tad5 = 1.248;   free.f_tad5 = 1;  units.f_tad5 = '-';       label.f_tad5 = 'scaled functional response for tadpole data set 5'; 
% par.f_tad6 = 1.725;   free.f_tad6 = 1;  units.f_tad6 = '-';       label.f_tad6 = 'scaled functional response for tadpole data set 6'; 
par.f_tad2 = 0.5896;    free.f_tad2 = 1;  units.f_tad2 = '-';       label.f_tad2 = 'scaled functional response for tadpole length weight data in AMA EPA tests';
par.f_tad3 = 1.253;    free.f_tad3 = 1;  units.f_tad3 = '-';       label.f_tad3 = 'scaled functional response for tadpole length weight data in AMAExt tests';
par.f_tad7 = 1.5;    free.f_tad7 = 1;  units.f_tad7 = '-';       label.f_tad7 = 'scaled functional response for tadpole length weight data in AMA test non public'; 

par.z_m = 2.651;       free.z_m   = 1;   units.z_m = '-';          label.z_m = 'zoom factor for males'; 
par.bn = 59.93;        free.bn = 1; units.bn='-'; label.bn='temperature stress scale';
%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
