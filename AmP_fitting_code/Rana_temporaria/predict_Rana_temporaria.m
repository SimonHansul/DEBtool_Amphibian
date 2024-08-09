function [prdData, info] = predict_Rana_temporaria(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  par.kV = 0.80;       % conversion efficiency from structure to reserve
  par.terminate=0;

  % compute temperature correction factors
  TC_ab = tempcorr(temp.ab, T_ref, [T_A, TrefH, TA_h]);
  TC = tempcorr(C2K(18), T_ref, [T_A, TrefH, TA_h]);
  TC22 = tempcorr(C2K(22), T_ref, [T_A, TrefH, TA_h]);
  TC19 = tempcorr(C2K(19), T_ref, [T_A, TrefH, TA_h]);
  
% zero-variate data

  % life cycle
  pars_tj = [g; k; l_T; v_Hb; v_Hj];               % compose parameter vector
  [t_j, t_b, l_j, l_b] = get_tp(pars_tj, f);       % -, scaled times & lengths at f
  pars_tp = [g; k; l_T; v_Hb; v_Hp];               % compose parameter vector
  [t_p, t_b, l_p, l_b, info] = get_tp(pars_tp, f); % -, scaled times & lengths at f
  
  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  % use the wet weight rescalinf 0.06 - 0.13 relation
  Ww_b = L_b^3 * d_V *(1 + f * w) / 0.06;      % g, wet weight at birth at f
  aT_b = t_b/ k_M/ TC_ab;           % d, age at birth at f and T

  % metam
  
  % puberty 
  tT_p = (t_p - t_j)/ k_M/ TC;      % d, time since metam at puberty at f and T

  % ultimate
  l_i = f - l_T;                    % -, scaled ultimate length
  L_i = L_m * l_i;                  % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;                % cm, ultimate physical length at f
  Ww_i = L_i^3 * (1 + f * w);       % g, ultimate wet weight 

  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
  RT_i = TC * reprod_rate(L_i, f, pars_R);                % #/d, ultimate reproduction rate at T

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC;                  % d, mean life span at T

  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve  
  E_0 = U_E0 * p_Am;     % J, initial reserve


  par.delM1 = del_Mt;
  par.delM2  = del_M;
  par.dV    = d_V;        % dry weight density (g/cm3)
  par.TA   = T_A;         % Arrhenius temperature (K)
  par.spAm = p_Am;        % max. surface-specific assimilation rate (J/(cm2 d))
  par.spM  = p_M;         % volume-specific somatic maintenance costs (J/(cm3 d))
  par.spT  = 0;           % surface-specific maintenance costs (J/(cm2 d))
  par.kJ   = k_J;         % maturity maintenance rate constant (1/d)
  par.EG   = E_G;         % volume-specific costs for growth (J/cm3)
  par.EHb  = E_Hb;        % maturity level at birth (J)
  par.EHj  = E_Hj;        % maturity level at metamorphosis (J)
  par.EH42 = E_H42;       % maturity at start of the metamorph. climax
  par.EHh  = E_Hh;        % maturity at hatching
  par.EHp  = E_Hp;        % maturity level at puberty (J)
  par.v    = v;           % energy conductance (cm/d)
  par.kap  = kap;         % allocation fraction to soma (-)
  par.kap2 = kap;         % minimum of the kappa value (to be discussed if it is useful or not)
  par.f2 = 0;             % starvation during climax
  par.kapR = kap_R ;      % reproduction efficiency (-)
  par.f_mo    = 1;        % scaled food density (-) from mother as initial condition
  par.yP   = 0.64;        % used only for starvation case, efficiency from reaserve to structure and from structure to reserve
  par.f = f;              % scaled functional response

  % temperature related parameters
  par.Texp = T_ref; % initial calculation to be done at reference temperature then we can retransform
  par.Tref = T_ref; % reference temperature of AmP values
  t_tmp = [0:0.1:300];
  X0(1,:) = E_0;     % state 1 is the reserve
  X0(2,:) = 0;      % state 3 maturity
  X0(3,:) = 1e-5;   % state 4 is the structural length
  X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
  X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
  X0(6,:) = 1;      % state 7 is survival probability (NOT USED)

  [to, Xo, par, phys,TE,YE, IE] = solve_deri(t_tmp, X0, par);
  aT_h2 = TE(1)/ TC_ab; % hatch time
  aT_b2 = TE(2)/ TC_ab;

  duration = (TE(4) - TE(2))/TC; % time from G25 to G46


  % pack to output
  prdData.ab = aT_b2;
  prdData.tj = duration; 
  prdData.tp = tT_p;
  prdData.am = aT_m;
  prdData.Lj = YE(4,3)/del_M;
  prdData.Li = Lw_i;
  prdData.Wwb = Ww_b;
  prdData.Wwi = Ww_i;
  prdData.Ri = RT_i;
  
  % univariate data
  % length-change in length for frogs
  rT_B = TC * k_M/ 3/ (1 + f_LdL/ g); L_i = f_LdL * L_m; 
  EdLw = rT_B * (L_i/del_Msvl - data.LdL(:,1));
  
  prdData.LdL = EdLw;


%% Ruthsatz data
  par.Texp = C2K(22); % initial calculation to be done at reference temperature then we can retransform
  par.Tref = T_ref; % reference temperature of AmP values
  pT_Am = TC22 * p_Am; vT = TC22 * v; pT_M = TC22 * p_M; kT_J = TC22 * k_J; 
  %E_0 = J_E_Am * U_E0 * mu_E; % J, initial reserve for kap constant
  par.f = fRu;
  par.spAm = pT_Am; 
  par.v = vT; 
  par.spM = pT_M ; 
  par.kJ = kT_J;
  %par.terminate=1;
  t_tmp = [0:0.1:300];
  X0(1,:) = E_0;     % state 1 is the reserve
  X0(2,:) = 0;      % state 3 maturity
  X0(3,:) = 1e-5;   % state 4 is the structural length
  X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
  X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
  X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
  X0_in = X0;

  hatchtime=aT_h2*TC_ab/TC22;
  [to, Xo, par, phys,TE,YE, IE] = solve_deri([0; data.massd(:,1)+hatchtime], X0, par);
  tha = TE(1);

  dw = (Xo(2:end,3).^3) .* d_V .* (1 + w .* phys.e(2:end) );
  prdData.massd = dw ./ arrayfun(@(x) rescaleweight(x, par.EH42, par.EHj), Xo(2:end,2));
  
  prdData.SVLd = (Xo(2:end,3)) ./ arrayfun(@(x) shape_coeff(x, par.EH42, par.EHj, del_Msvl, par.delM2), (Xo(2:end,2)));    %del_Msvl;

  prdData.dstage = [TE(2:4)-tha]; % duration of the stages from hatch

  tempdata = (TE(3) - tha) * TC22;
  prdData.dtemp = (tempdata ./ tempcorr(C2K(data.dtemp(:,1)),T_ref, [T_A, TrefH, TA_h])); 
  
  svlstemp22 = 10. * YE(3,3)/del_Mt;
  masstemp22 = YE(3,3)^3 * d_V * (1 + w * fRu ) ./ 0.06; % wet weight at the beginning of the metamorphosis

  % additional datasets to increase the weight of the time spots
  prdData.tbirth=TE(2)-TE(1);
  prdData.tmetsta=TE(3)-TE(1);
  prdData.tmetsto=TE(4)-TE(1);

%% Ruthsatz temperature dependent data for weight and lengths
svlstemp = zeros(length(data.dsvlltemp(:,1)),1);
totltemp = zeros(length(data.dtotltemp(:,1)),1);
masstemp = zeros(length(data.dtotltemp(:,1)),1);
durtemp  = zeros(length(data.dtotltemp(:,1)),1);
for i=1:length(data.dtotltemp(:,1))
  par.f = fRu;
  par.bn=bn;
  par.Texp = C2K(data.dtotltemp(i,1)); % initial calculation to be done at reference temperature then we can retransform
  par.Tref = T_ref; % reference temperature of AmP values
  TCtemp=tempcorr(par.Texp, T_ref, [T_A, TrefH, TA_h]);
  pT_Am = TCtemp * p_Am; vT = TCtemp * v; pT_M = TCtemp * p_M; kT_J = TCtemp * k_J; 
  par.spAm = pT_Am; 
  par.v = vT; 
  par.spM = pT_M ; 
  par.kJ = kT_J;
  t_tmp = [0:0.1:500]';
  X0(1,:) = E_0;     % state 1 is the reserve
  X0(2,:) = 0;      % state 3 maturity
  X0(3,:) = 1e-5;   % state 4 is the structural length
  X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
  X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
  X0(6,:) = 1;      % state 7 is survival probability (NOT USED)

  hatchtime=aT_h2*TC_ab/TCtemp;
  par.terminate=2; % stop when reaching start of the metamorphosis climax
  [to, Xo, par, phys,TE,YE,IE] = solve_deri(t_tmp, X0, par);
  durtemp(i)  = TE(3) - TE(1);
  svlstemp(i) =  YE(3,3)/del_Msvl;
  totltemp(i) =  YE(3,3)/del_Mt;
  masstemp(i) = YE(3,3)^3 * d_V * (1 + w .* fRu ) / 0.06; % wet weight at the beginning of met. climax
end 
prdData.dtemp = durtemp;
prdData.dsvlltemp = svlstemp;
prdData.dtotltemp = totltemp;
prdData.dmasstemp = masstemp;

%% Morand97 data
masstemp = zeros(length(data.tempsdur(:,1)),1);
durtemp  = zeros(length(data.tempsdur(:,1)),1);
for i=1:length(data.tempsdur(:,1))
  par.f = f;
  par.bn=bn;
  par.Texp = C2K(data.tempsdur(i,1)); % initial calculation to be done at reference temperature then we can retransform
  par.Tref = T_ref; % reference temperature of AmP values
  TCtemp=tempcorr(par.Texp, T_ref, T_A);
  pT_Am = TCtemp * p_Am; vT = TCtemp * v; pT_M = TCtemp * p_M; kT_J = TCtemp * k_J; 
  par.spAm = pT_Am; 
  par.v = vT; 
  par.spM = pT_M ; 
  par.kJ = kT_J;
  t_tmp = [0:0.1:500]';
  X0(1,:) = E_0;     % state 1 is the reserve
  X0(2,:) = 0;      % state 3 maturity
  X0(3,:) = 1e-5;   % state 4 is the structural length
  X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
  X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
  X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
  par.terminate=1;
  [to, Xo, par, phys,TE,YE, IE] = solve_deri(t_tmp, X0, par);
  durtemp(i)  = TE(4)-TE(1);
  masstemp(i) = YE(4,3)^3 * d_V * (1 + w .* f ); % wet weight at the beginning of met. climax
end 
prdData.tempsdur = durtemp;
prdData.tempdry = masstemp;

  %% Aviles22 data
  % here a linear relation was assumed between stage and time to reach
  % metamorphosis
  TC16p5 = tempcorr(C2K(16.5), T_ref, [T_A, TrefH, TA_h]);
  par.Texp = C2K(16.5);
  pT_Am = TC16p5 * p_Am; vT = TC16p5 * v; pT_M = TC16p5 * p_M; kT_J = TC16p5 * k_J; 
  par.spAm = pT_Am; 
  par.v = vT; 
  par.spM = pT_M ; 
  par.kJ = kT_J;
  par.f = fAv;
  X0(1,:) = E_0;     % state 1 is the reserve
  X0(2,:) = 0;      % state 3 maturity
  X0(3,:) = 1e-5;   % state 4 is the structural length
  X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
  X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
  X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
  par.terminate=1;
  birthtime=aT_b2*TC_ab/TC16p5;
  meanstage=42.4;
  if birthtime > 48
      prdData.svltad = 0. * data.svltad(:,2);
      prdData.totalentad = 0. * data.totalentad(:,2);
      prdData.weighttad = 0. * data.weighttad(:,2);
      prdData.startj = 0;
      return
  end

  timeaxis = [0:0.1:100];
  [to, Xo, par, phys,TE,YE,IE] = solve_deri(timeaxis, X0, par);
  t42=TE(3);
  prdData.startj = t42;
  prdData.ww42 = (  YE(3,3).^3 .* d_V .* (1 + w .* fAv)) / 0.06;
  prdData.svl42 =  YE(3,3) / par.del_Msvl;
  prdData.totlen42 =  YE(3,3) / par.delM1;






  
  
  


  