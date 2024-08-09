function [prdData, info] = predict_Bufo_bufo(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  par.kV = 0.80;       % conversion efficiency from structure to reserve
  par.terminate=0;



  % compute temperature correction factors
  TC_ab = tempcorr(temp.ab, T_ref, T_A);
  TC    = tempcorr(temp.am, T_ref, T_A);
  TC20p5 = tempcorr(C2K(temp.tL), T_ref, T_A);
  TC22 = tempcorr(C2K(22), T_ref, T_A);
  TC19 = tempcorr(C2K(19), T_ref, T_A);
  TC18 = tempcorr(C2K(18), T_ref, T_A);
  % zero-variate data

%% solution with the amphibian model
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve  
  E_0 = U_E0 * p_Am;     % J, initial reserve


  par.delM1 = del_Msvl;  % SVL shape coefficient
  par.delM2  = del_M;    % shape coefficient toad (exppected to be very similar to svl)
  par.dV    = d_V;       % dry weight density (g/cm3)
  par.TA   = T_A;        % Arrhenius temperature (K)
  par.spAm = p_Am;       % max. surface-specific assimilation rate (J/(cm2 d))
  par.spM  = p_M;        % volume-specific somatic maintenance costs (J/(cm3 d))
  par.spT  = 0;          % surface-specific maintenance costs (J/(cm2 d))
  par.kJ   = k_J;        % maturity maintenance rate constant (1/d)
  par.EG   = E_G;        % volume-specific costs for growth (J/cm3)
  par.EHb  = E_Hb;       % maturity level at birth (J)
  par.EHj  = E_Hj;       % maturity level at metamorphosis (J)
  par.EH42 = E_H42;      % maturity at start of the metamorph. climax
  par.EHh  = E_Hh;       % maturity at hatching
  par.EHp  = E_Hp;       % maturity level at puberty (J)
  par.v    = v;          % energy conductance (cm/d)
  par.kap  = kap;        % allocation fraction to soma (-)
  par.kap2 = kap;        % minimum of the kappa value (to be discussed if it is useful or not)
  par.f2 = 0;            % starvation during climax
  par.kapR = kap_R ;     % reproduction efficiency (-)
  par.f_mo    = 1;       % scaled food density (-) from mother as initial condition
  par.yP   = 0.64;       % used only for starvation case, efficiency from reaserve to structure and from structure to reserve
  par.f = f;             % scaled functional response

  % temperature related parameters
  par.Texp = temp.ab; % initial calculation to be done at the typical temperature of 15 degrees
  par.Tref = T_ref; % reference temperature of AmP values
  
  t_tmp = [0:0.1:1000];
  X0(1,:) = E_0;     % state 1 is the reserve
  X0(2,:) = 0;       % state 3 maturity
  X0(3,:) = 1e-5;    % state 4 is the structural length
  X0(4,:) = 0;       % state 5 is the cumulative reproduction (NOT USED)
  X0(5,:) = 0;       % state 6 is the scaled damage (referenced to water) (NOT USED)
  X0(6,:) = 1;       % state 7 is survival probability (NOT USED)

  [to, Xo, parout, phys,TE,YE, IE] = solve_deri(t_tmp, X0, par);
  aT_b2 = TE(2);  % age at birth;

  % life cycle 
  pars_tj = [g k l_T v_Hb v_Hj]; 
  [t_j, t_b, l_j, l_b] = get_tp(pars_tj, f);
  pars_tp = [g; k; l_T; v_Hb; v_Hp]; % compose parameter vector
  [t_p, t_b, l_p, l_b, info] = get_tp(pars_tp, f); % -, scaled length at birth at f

  % birth
  L_b = L_m * l_b;                % cm, structural length at birth at f
  aT_b = t_b/ k_M/ TC_ab;            % d, age at birth at f and T
  
  % puberty
  L_p = L_m * l_p;               % cm, structural length at puberty at f
  Lw_p = L_p/ del_M;             % cm, body length at puberty at f
  tT_p = (t_p - t_j)/ k_M/ TC;   % d, time time metam at puberty at f and T
  
  % ultimate for toad
  l_i = f - l_T;                 % -, scaled ultimate length at f
  L_i = L_m * l_i;               % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;             % cm, ultimate body length at f
  Ww_i = L_i^3 * (1 + f * w);    % g, ultimate wet weight 
  
  % reproduction
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
  RT_i = TC * reprod_rate(L_i, f, pars_R);                  % #/d, ultimate reproduction rate at T

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G]; % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);     % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC;                 % d, mean life span at T
  
  % males
  p_Am_m = z_m * p_M/ kap;             % J/d.cm^2, {p_Am} spec assimilation flux
  E_m_m = p_Am_m/ v;                   % J/cm^3, reserve capacity [E_m]
  g_m = E_G/ (kap* E_m_m);             % -, energy investment ratio
  L_mm = v/ k_M/ g_m;                  % cm, max struct length
  m_Em_m = y_E_V * E_m_m/ E_G;         % mol/mol, reserve capacity 
  w_m = m_Em_m * w_E/ w_V;             % -, contribution of reserve to weight
  pars_tpm = [g_m k l_T v_Hb v_Hp];
  [t_pm, t_bm, l_pm, l_bm] = get_tp(pars_tpm, f);
  L_im = f * L_mm; Lw_im = L_im/ del_M;    % cm, ultimate SVL
  L_pm = l_pm * L_mm; Lw_pm = L_pm/ del_M; % cm,  SVL at puberty
  Ww_im = L_im^3 * (1 + f * w_m); % g, ultimate weight



  % pack to output
  prdData.ab = aT_b2;
  
  prdData.tp = tT_p;
  prdData.am = aT_m;
  prdData.Lj = YE(4,3)/del_M;
  prdData.Lp = Lw_p;
  prdData.Li = Lw_i;
  prdData.Lim = Lw_im;
  prdData.Wwj = YE(4,3)^3 * d_V * (1 + f * w) / 0.13;
  prdData.Wwi = Ww_i;
  prdData.Wwim = Ww_im;
  prdData.Ri = RT_i;
  
  % uni-variate data
  %% katzmann data
  par.Texp = C2K(20.5); % initial calculation to be done at reference temperature then we can retransform
  par.Tref = T_ref; % reference temperature of AmP values
  pT_Am = TC20p5 * p_Am; vT = TC20p5 * v; pT_M = TC20p5 * p_M; kT_J = TC20p5 * k_J; 
  %E_0 = J_E_Am * U_E0 * mu_E; % J, initial reserve for kap constant
  par.f=f_tL;
  par.spAm = pT_Am; 
  par.v = vT; 
  par.spM = pT_M ; 
  par.kJ = kT_J;

  t_tmp = [0:0.1:600];
  X0(1,:) = E_0;     % state 1 is the reserve
  X0(2,:) = 0;      % state 3 maturity
  X0(3,:) = 1e-5;   % state 4 is the structural length
  X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
  X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
  X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
  X0_in = X0;
  par.terminate=1;
  [to, Xo, par, phys,TE,YE, IE] = solve_deri(t_tmp, X0, par);
  Xb = YE(2,:);
  X42 = YE(3,:);
  Xj = YE(4,:);
  
  prdData.tj = TE(4); % event of finishing metamorphosis

  par.terminate =0;
  [to, Xo, par, phys,TE,YE, IE] = solve_deri([0; data.tL(:,1) ], Xb, par);
  Lvals = Xo(2:end,3);
  EHvals = Xo(2:end,2);
  if length(to)<=length(data.tL(:,1))
      Lvals = [Lvals; zeros(length(data.tL(:,1))-length(to)+1,1)];
      EHvals = [Xo(2:end,2); zeros(length(data.tL(:,1))-length(to)+1,1)];
  end
  prdData.tL = Lvals ./ arrayfun(@(x) shape_coeff(x, par.EH42, par.EHj, del_Msvl, par.delM2), EHvals);    %del_Msvl;

  prdData.SL = 10*[Xb(3)/del_Msvl; X42(3)/del_Msvl; Xj(3)/del_M];

  prdData.SW = [X42(3)^3 * d_V * (1 + f_tL * w) / 0.06; Xj(3)^3 * d_V * (1 + f_tL * w) / 0.13];

  
%% Brunellli data
  par.Texp = C2K(22); % initial calculation to be done at reference temperature then we can retransform
  par.Tref = T_ref; % reference temperature of AmP values
  pT_Am = TC22 * p_Am; vT = TC22 * v; pT_M = TC22 * p_M; kT_J = TC22 * k_J; 
  %E_0 = J_E_Am * U_E0 * mu_E; % J, initial reserve for kap constant
  par.f= fBr;
  par.spAm = pT_Am; 
  par.v = vT; 
  par.spM = pT_M ; 
  par.kJ = kT_J;
  t_tmp = [0:0.1:600];
  X0(1,:) = E_0;     % state 1 is the reserve
  X0(2,:) = 0;      % state 3 maturity
  X0(3,:) = 1e-5;   % state 4 is the structural length
  X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
  X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
  X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
  par.termiante = 1;
  [~,~, ~, ~,TE,~, ~] = solve_deri(t_tmp, X0, par);
  par.terminate =0;
  [~, Xo1, par, ~,~,~, ~] = solve_deri([0; data.tW(:,1) + TE(2)], X0, par);
  dw = (Xo1(2:end,3)).^3 * d_V * (1 + fBr * w);
  prdData.tW = dw ./ arrayfun(@(x) rescaleweight(x, par.EH42, par.EHj), Xo1(2:end,2));
  prdData.St = [TE(3)-TE(2); TE(4)-TE(2)];

%% Morand97 data
  par.Texp = C2K(15); % initial calculation to be done at reference temperature then we can retransform
  par.Tref = T_ref; % reference temperature of AmP values
  TC15 = tempcorr(C2K(15),T_ref, T_A);
  TC27 = tempcorr(C2K(27),T_ref, T_A);
  pT_Am = TC15 * p_Am; vT = TC15 * v; pT_M = TC15 * p_M; kT_J = TC15 * k_J; 
  %E_0 = J_E_Am * U_E0 * mu_E; % J, initial reserve for kap constant
  par.f = f;
  par.spAm = pT_Am; 
  par.v = vT; 
  par.spM = pT_M ; 
  par.kJ = kT_J;

  X0(1,:) = E_0;     % state 1 is the reserve
  X0(2,:) = 0;      % state 3 maturity
  X0(3,:) = 1e-5;   % state 4 is the structural length
  X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
  X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
  X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
  X0_in = X0;
%   par.kap=0.95;
%   par.EHj = 95;
  %par.kV = 0.99;
  par.terminate=1;
  [to, Xo, par, phys,TE,YE, IE] = solve_deri(t_tmp, X0, par);
  dur46_15 = TE(4)-TE(2);
  Lv=YE(4,3);
  E_M = pT_Am / vT;
  evl = YE(4,1) ./ (E_M * Lv.^3) ;
  dryw46_15 = 1000. * Lv^3 * d_V * (1+w*evl); % factor 1000 because data in mg

  par.Texp = C2K(27);
  pT_Am = TC27 * p_Am; vT = TC27 * v; pT_M = TC27 * p_M; kT_J = TC27 * k_J; 
  par.spAm = pT_Am; 
  par.v = vT; 
  par.spM = pT_M ; 
  par.kJ = kT_J;

  X0(1,:) = E_0;     % state 1 is the reserve
  X0(2,:) = 0;      % state 3 maturity
  X0(3,:) = 1e-5;   % state 4 is the structural length
  X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
  X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
  X0(6,:) = 1;      % state 7 is survival probability (NOT USED)

  par.terminate=1; % stop at end of metamorphosis
  [to, Xo, par, phys,TE,YE,IE] = solve_deri(t_tmp, X0, par);
  dur46_27 = TE(4)-TE(2);
  Lv=YE(4,3);
  E_M = pT_Am / vT;
  evl = YE(4,1) ./ (E_M * Lv.^3) ;
  dryw46_27 = 1000. * Lv^3 * d_V * (1+w*evl); 

  prdData.tempsdur=[dur46_15;dur46_27];
  prdData.tempdry =[dryw46_15;dryw46_27];

%% Laurila data
par.Texp = C2K(19); % initial calculation to be done at reference temperature then we can retransform
par.Tref = T_ref; % reference temperature of AmP values
pT_Am = TC19 * p_Am; vT = TC19 * v; pT_M = TC19 * p_M; kT_J = TC19 * k_J; 
%E_0 = J_E_Am * U_E0 * mu_E; % J, initial reserve for kap constant

%low food
par.f = fL1;
par.spAm = pT_Am; 
par.v = vT; 
par.spM = pT_M ; 
par.kJ = kT_J;
t_tmp = [0:0.1:300];
X0(1,:) = E_0;     % state 1 is the reserve
X0(2,:) = 0;      % state 3 maturity
X0(3,:) = 1e-5;   % state 4 is the structural length
X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
X0_in = X0;
par.terminate=1;
[to, Xo, par, phys,TE,YE, IE] = solve_deri(t_tmp, X0, par);


prdData.fD = TE(3)-TE(2);
prdData.fW = 1000. * YE(3,3)^3 * d_V * (1+fL1*w) / 0.06;

% high food
par.f = fL2;
par.spAm = pT_Am; 
par.v = vT; 
par.spM = pT_M ; 
par.kJ = kT_J;
t_tmp = [0:0.1:300];
X0(1,:) = E_0;     % state 1 is the reserve
X0(2,:) = 0;      % state 3 maturity
X0(3,:) = 1e-5;   % state 4 is the structural length
X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
X0_in = X0;
par.terminate=1;
[to, Xo, par, phys,TE,YE, IE] = solve_deri(t_tmp, X0, par);


prdData.fD2 = TE(3)-TE(2);
prdData.fW2 = 1000. * YE(3,3)^3 * d_V * (1+fL2*w) / 0.06;


%% Mikó data
par.Texp = C2K(18); % initial calculation to be done at reference temperature then we can retransform
par.Tref = T_ref; % reference temperature of AmP values
pT_Am = TC18 * p_Am; vT = TC18 * v; pT_M = TC18 * p_M; kT_J = TC18 * k_J; 

%low food
par.f = 1;
par.spAm = pT_Am; 
par.v = vT; 
par.spM = pT_M ; 
par.kJ = kT_J;
t_tmp = [0:0.1:300];
X0(1,:) = E_0;     % state 1 is the reserve
X0(2,:) = 0;      % state 3 maturity
X0(3,:) = 1e-5;   % state 4 is the structural length
X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
X0_in = X0;
par.terminate=1;
[to, Xo, par, phys,TE,YE, IE] = solve_deri(t_tmp, X0, par);


prdData.tmet = TE(3)-TE(2);
prdData.Mmet = 1000. * YE(3,3)^3 * (1+fL1*w) * d_V / 0.06;
end

