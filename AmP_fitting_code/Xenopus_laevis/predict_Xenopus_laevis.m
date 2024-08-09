function [prdData, info] = predict_Xenopus_laevis(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

  % % customized filters to contrain additional parameter

  %filterChecks =  fWa > 1.5 || fWa < 0.5 || fNa > 1.1 || fNa < 0.5 || f_tad7 > 1.5 || f_tad7 < 0.5 || shiftNa <= 0;
  filterChecks =  fNa > 1.1 || fNa < 0.5 || f_tad7 > 1.5 || f_tad7 < 0.5 || shiftNa <= 0 || TrefH <= T_ref*1.01;

  if filterChecks
      info = 0;
      prdData = {};
      return;
  end
  
  %% pass the parameter to a new par object
  parsolver.delM1 = del_Mt;
  parsolver.delM2  = del_M;
  parsolver.d_V    = d_V;        % dry weight density (g/cm3)
  parsolver.TA   = T_A;         % Arrhenius temperature (K)
  parsolver.spAm = p_Am;     % max. surface-specific assimilation rate (J/(cm2 d))
  parsolver.spM  = p_M;     % volume-specific somatic maintenance costs (J/(cm3 d))
  parsolver.spT  = 0;            % surface-specific maintenance costs (J/(cm2 d))
  parsolver.kJ   = k_J;        % maturity maintenance rate constant (1/d)
  parsolver.EG   = E_G;       % volume-specific costs for growth (J/cm3)
  parsolver.EHb  = E_Hb;       % maturity level at birth (J)
  parsolver.EHj  = E_Hj;      % maturity level at metamorphosis (J)
  parsolver.EH42 = E_H62;     % maturity at start of the metamorph. climax
  parsolver.EHh  = E_Hh;       % maturity at hatching
  parsolver.EHp  = E_Hp;   % maturity level at puberty (J)
  parsolver.v    = v;       % energy conductance (cm/d)
  parsolver.kap  = kap;         % allocation fraction to soma (-)
  parsolver.kap2 = kap;         % minimum of the kappa value (to be discussed if it is useful or not)
  parsolver.f2 = 0;              % starvation during climax
  parsolver.kapR = kap_R ;        % reproduction efficiency (-)
  parsolver.f_mo    = 1;         % scaled food density (-) from mother as initial condition
  parsolver.yP   = 0.64;
  parsolver.f = f;
  parsolver.bn = bn;
  parsolver.dV = d_V;
  parsolver.terminate = 0;
  parsolver.kV = 0.8;
 


  % compute temperature correction factors
  TC = tempcorr(temp.ab, T_ref, [T_A TrefL TrefH TA_l TA_h]);
  TC_am = tempcorr(temp.am, T_ref, [T_A TrefL TrefH TA_l TA_h]);
  TC_Ri = tempcorr(temp.Ri, T_ref, [T_A TrefL TrefH TA_l TA_h]);
  TC_tL = tempcorr(temp.tL, T_ref, [T_A TrefL TrefH TA_l TA_h]);
  TC_tLtad = tempcorr(temp.tLtad, T_ref, [T_A TrefL TrefH TA_l TA_h]);
  TC_tAMA = tempcorr(C2K(22), T_ref, [T_A TrefL TrefH TA_l TA_h]);
  TC_tLmet = tempcorr(temp.tLmet, T_ref, [T_A TrefL TrefH TA_l TA_h]);
  

  % hatch
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve  
  E_0 = U_E0 * p_Am;     % J, initial reserve
  Wd0 = E_0 * w_E/ mu_E ; % g, egg dry weight 
  [U_H, aUL] = ode45(@dget_aul, [0; U_Hh; U_Hb], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
 
  % zero-variate date
  parsolver.Texp = C2K(23); % initial calculation to be done at reference temperature then we can retransform
  parsolver.Tref = T_ref; % reference temperature of AmP values
  pT_Am = TC * p_Am;
  vT = TC * v;
  pT_M = TC * p_M; 
  kT_J = TC * k_J; 
  parsolver.spAm = pT_Am; 
  parsolver.v = vT; 
  parsolver.spM = pT_M ; 
  parsolver.kJ = kT_J;
  parsolver.f = f_tL;
  t_tmp = [0:0.1:1000];
  X0(1,:) = E_0;     % state 1 is the reserve
  X0(2,:) = 0;      % state 3 maturity
  X0(3,:) = 1e-5;   % state 4 is the structural length
  X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
  X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
  X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
  X0_in = X0;
  [to, Xo, parout, phys,TE,YE] = solve_deri(t_tmp, X0, parsolver);
  aT_h = TE(1);%/ TC;
  aT_b = TE(2);%/ TC;
  aT_62 = TE(3);%/ TC;
  if (length(TE)>3)
     aT_j = TE(4);%/ TC;
     if (length(TE)>4)
        aT_p = TE(5);%/ TC; 
     else
        aT_p = 2000;
     end
  else
      aT_j = 1000;
      aT_p = 2000;
  end
  tT_j = aT_j - aT_b;
  tT_p = aT_p - aT_j;
  
  % life cycle for tadpole
  pars_tpj = [g k l_T v_Hb v_Hj]; % for tadpole
  [t_j, t_b, l_j, l_b, info] = get_tp(pars_tpj, f);
  pars_tp = [g k l_T v_Hb v_Hp]; % for frog
  [t_p, t_b, l_p, l_b] = get_tp(pars_tp, f);
    

  % birth
  L_b = L_m * l_b;                  % cm, structural length at birth at f
  Ww_b = L_b^3 * d_V * (1 + f * w) / 0.06;       % g, wet weight at birth at f. Corrected according to weight rescaling
  
  % ultimate for toad
  l_i = f - l_T;                  % -, scaled ultimate length at f
  L_i = L_m * l_i;                % cm, ultimate structural length at f
  Lw_i = L_i/ del_M;              % cm, ultimate total length at f
  Ww_i = L_i^3 * (1 + f * w);     % g, ultimate wet weight 
  
  % reproduction
  % does not need to be recomputed with the new code, because the
  % differences with the solution without the metamorphosis are negligible
  % when dealing with the reproduction rate
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
  RT_i = TC_Ri * reprod_rate(L_i, f, pars_R);               % #/d, ultimate reproduction rate at T

  % life span
  pars_tm = [g; l_T; h_a/ k_M^2; s_G]; % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);     % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;              % d, mean life span at T

  % males
  p_Am_m = z_m * p_M/ kap;             % J/d.cm^2, {p_Am} spec assimilation flux
  E_m_m = p_Am_m/ v;                   % J/cm^3, reserve capacity [E_m]
  g_m = E_G/ (kap* E_m_m);             % -, energy investment ratio
  m_Em_m = y_E_V * E_m_m/ E_G;         % mol/mol, reserve capacity 
  w_m = m_Em_m * w_E/ w_V;             % -, contribution of reserve to weight
  L_mm = v/ k_M/ g_m;                  % cm, max struct length
  L_im = f * L_mm; Lw_im = L_im/ del_M;% cm, ultimate structural, body length
  Ww_im = L_im^3 * (1 + f * w_m);      % g, ultimate wet weight

  % pack to output
  prdData.ah = aT_h;
  prdData.ab = aT_b;
  prdData.tj = tT_j;
  prdData.tp = tT_p;
  prdData.am = aT_m;
  prdData.Li = Lw_i;
  prdData.Lim = Lw_im;
  prdData.Wwb = Ww_b;
  prdData.Wwi = Ww_i;
  prdData.Wwim = Ww_im;
  prdData.Ri = RT_i;
  
  %%% uni-variate data

  %% time - length for frogs (past metamorphosis)
  parsolver.Texp = C2K(20);
  X0 = X0_in;
  pT_Am = TC_tL * p_Am;
  vT = TC_tL * v;
  pT_M = TC_tL * p_M; 
  kT_J = TC_tL * k_J; 
  E_0 = J_E_Am * U_E0 * mu_E; % J, initial reserve for kap constant
  X0(1) = E_0;
  parsolver.spAm = pT_Am; 
  parsolver.v = vT; 
  parsolver.spM = pT_M ; 
  parsolver.kJ = kT_J;
  parsolver.f = f_tL;
  parsolver.terminate=0;
  [to, Xo, parout, phys,TE,YE] = solve_deri([0; tL(:,1)+( aT_j * TC/TC_tL)], X0, parsolver);
  % t-L data for frogs
  rT_B = TC_tL * k_M/ 3/ (1 + f_tL/ g); % d 1/von Bert growth rate
  L_i = f_tL * L_m;  
  L_j = YE(4,3);
%   ELw = Xo(2:end,3)/del_M;
%   EWw = Xo(2:end,3).^3 * (1 + f_tL * w);
  ELw2 = (L_i - (L_i - L_j) * exp( - rT_B * (tL(:,1) )))/ del_M; % cm, body length at time
  EWw2 = (L_i - (L_i - L_j) * exp( - rT_B * (tWw(:,1) ))).^3 * (1 + f_tL * w); % g, wet weight  

  % pack to output
  prdData.tL = ELw2;
  prdData.tWw = EWw2;



  %% data for the tadpoles (NF book)
  parsolver.terminate=0; % solve the ODE even after the end of the metamorphosis
  parsolver.Texp = C2K(23);
  pT_Am = TC_tLtad * p_Am; 
  vT = TC_tLtad * v; 
  pT_M = TC_tLtad * p_M; 
  kT_J = TC_tLtad * k_J; 
  parsolver.spAm = pT_Am; 
  parsolver.v = vT; 
  parsolver.spM = pT_M ; 
  parsolver.kJ = kT_J;
  parsolver.f = f_tad;
  [t, Xout, parsolver, phys,~,~] = solve_deri(tLtad(:,1), X0, parsolver);
  ELw_tad = Xout(:,3) ./ arrayfun(@(x) shape_coeff(x, parsolver.EH42, parsolver.EHj, parsolver.delM1, parsolver.delM2), Xout(:,2)); % physical length back in mm
  parsolver.terminate=1;
  [t, Xout, parsolver, phys,TE,YE] = solve_deri(t_tmp, X0, parsolver);
  prdData.tLtad1 = 10. .* YE(1:3,3) ./ arrayfun(@(x) shape_coeff(x, parsolver.EH42, parsolver.EHj, parsolver.delM1, parsolver.delM2), YE(1:3,2));
  prdData.tLtad2 = TE(1:3);  
  prdData.tLtad = ELw_tad;

%% Extended AMA time dependent quantities (3 different datasets)
  parsolver.terminate=2; % stop the ODE even at beginning of climax
  parsolver.Texp = C2K(22);
  pT_Am = TC_tAMA * p_Am; 
  vT = TC_tAMA * v; 
  pT_M = TC_tAMA * p_M; 
  kT_J = TC_tAMA * k_J; 
  %E_0 = J_E_Am * U_E0 * mu_E; % J, initial reserve for kap constant
  parsolver.spAm = pT_Am; 
  parsolver.v = vT; 
  parsolver.spM = pT_M ; 
  parsolver.kJ = kT_J;
  parsolver.f = f_tad7;

  [gc,gr] = groupcounts(data.tLtadExt(:,1));
  t_tmp = union(t_tmp, data.tLtadExt(:,1) + 9);  % 9 should be the time from spawn to start of test for this dataset
  [t, Xout, parsolver, phys,TE,YE] = solve_deri(t_tmp, X0, parsolver);
  idx = find(ismember(t, data.tLtadExt(:,1) + 9));  % 9 should be the time from spawn to start of test
  valL = Xout(idx,3) ./ arrayfun(@(x) shape_coeff(x, parsolver.EH42, parsolver.EHj, del_Me, parsolver.delM2), Xout(idx,2));
  valdw = Xout(idx,3).^3 .* d_V .* (1 + phys.e(idx) .* w); % compute first dry weight
  valW = valdw ./ arrayfun(@(x) rescaleweight(x, parsolver.EH42, parsolver.EHj), Xout(idx,2));
  valL2=[];
  valW2=[];
  for i=1:length(idx)
      valL2=[valL2; repmat(valL(i),gc(i), 1)];
      valW2=[valW2; repmat(valW(i),gc(i), 1)];
  end
  prdData.tLtadExt = valL2;
  prdData.tWtadExt = valW2;

  % SVL - weight
  EWw_tad = (LWtad(:,1) * del_Me).^3 * d_V * (1 + f_tad2 * w) ./ 0.06;  % g, wet weight, assumption that are all pre-climax FOR EPA data
  EWw_tad1 = (LWtad1(:,1) * del_Me).^3 * d_V * (1 + f_tad3 * w)./ 0.06;  % g, wet weight, assumption that are all pre-climax
  
  % Hind Limb Length - weight
  EWw_tad2 = (LWtad2(:,1) * del_Meh).^3 * d_V * (1 + f_tad3 * w)./ 0.06;  % g, wet weight, assumption that are all pre-climax
 
  % pack to output
  prdData.LWtad = EWw_tad;
  prdData.LWtad1 = EWw_tad1;
  prdData.LWtad2 = EWw_tad2;


  %% Temperature - time for emrbyos
  % this to be kept this way because it refers to embryonic stage, where
  % the temperature should not be inducing a change in the kappa value.
  ETt = 1 ./ tempcorr(C2K(data.Tt(:,1)), C2K(25), [T_A TrefL TrefH TA_l TA_h]); 
  prdData.Tt = ETt;


  %% Nations'11 data during met. development
  pT_Am = TC_tLmet * p_Am; vT = TC_tLmet * v; pT_M = TC_tLmet * p_M; kT_J = TC_tLmet * k_J; 
  E_0 = J_E_Am * U_E0 * mu_E; % J, initial reserve for kap constant
  X0(1,:) = E_0;
  X0(2,:) = 0;      % state 3 maturity
  X0(3,:) = 1e-5;   % state 4 is the structural length
  X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
  X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
  X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
  parsolver.f = fNa;
  parsolver.spAm = pT_Am; 
  parsolver.v = vT; 
  parsolver.spM = pT_M ; 
  parsolver.kJ = kT_J;
  parsolver.f = f;
  parsolver.Texp = C2K(23.6);
  parsolver.terminate=1;
  [t, Xout, ~, phys,TE,YE] = solve_deri(t_tmp, X0_in, parsolver);
  prdData.tLmet1 = 10. .* YE(3:4,3) ./ arrayfun(@(x) shape_coeff(x, parsolver.EH42, parsolver.EHj, parsolver.delM1, parsolver.delM2), YE(3:4,2));
  % shiftNa is fixed to a value of 2 in the pars_init file
  prdData.tLmet2 = (TE(3:4) - shiftNa);  
  parsolver.terminate = 0;
  [tout, Xout, parout, phys,~,~] = solve_deri([0; (tLmet(:,1) + shiftNa)], X0, parsolver);
  ELtmet = 10. * Xout(2:end,3) ./ arrayfun(@(x) shape_coeff(x, parsolver.EH42, parout.EHj, parsolver.delM1, parsolver.delM2), Xout(2:end,2));
  ELtmetsvl = 10. * Xout(2:end,3) ./ arrayfun(@(x) shape_coeff(x, parsolver.EH42, parout.EHj, del_Me, del_M ), Xout(2:end,2));
  prdData.tLmet = ELtmet;
  prdData.tLmetsvl = ELtmetsvl;


%% Ruthsatz data, on relative duration of pre-met. phase
    % calculation at 22 degrees
      parsolver.f = fRu;
      parsolver.bn=bn;
      parsolver.Texp = C2K(22); % initial calculation to be done at reference temperature then we can retransform
      parsolver.Tref = T_ref; % reference temperature of AmP values
      TCtemp=tempcorr(parsolver.Texp, T_ref, [T_A TrefL TrefH TA_l TA_h]);
      pT_Am = TCtemp * p_Am; vT = TCtemp * v; pT_M = TCtemp * p_M; kT_J = TCtemp * k_J; 
      parsolver.spAm = pT_Am; 
      parsolver.v = vT; 
      parsolver.spM = pT_M ; 
      parsolver.kJ = kT_J;
      t_tmp = [0:0.1:500];
      X0(1,:) = E_0;     % state 1 is the reserve
      X0(2,:) = 0;      % state 3 maturity
      X0(3,:) = 1e-5;   % state 4 is the structural length
      X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
      X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
      X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
      parsolver.terminate=2;
      [to, Xo, ~, phys,TE,YE,IE] = solve_deri(t_tmp, X0, parsolver);
      % ind42 = find(IE==3);
      durtemp22  = TE(3) - TE(2);
      svlstemp22 = 10. * YE(3,3)/del_Me;
      totltemp22 = 10. * YE(3,3)/del_Mt;
      masstemp22 = 1000.*YE(3,3)^3 * d_V * (1 + w * fRu )/ 0.06;

  %% Ruthsatz data on temperature stress
  svlstemp = zeros(length(data.dmetsvl(:,1)),1);
  totltemp = zeros(length(data.dmetsvl(:,1)),1);
  masstemp = zeros(length(data.dmetsvl(:,1)),1);
  durtemp  = zeros(length(data.dmetsvl(:,1)),1);
  X0(1,:) = E_0;     % state 1 is the reserve
  X0(2,:) = 0;      % state 3 maturity
  X0(3,:) = 1e-5;   % state 4 is the structural length
  X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
  X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
  X0(6,:) = 1;      % state 7 is survival probability (NOT USED)  
  temparr = data.dmetsvl(:,1);
  T_ref=par.T_ref;
  T_A=par.T_A;
  TrefL = par.TrefL;
  TrefH = par.TrefH;
  TA_l = par.TA_l;
  TA_h = par.TA_h;
  p_Am = p_Am; 
  v = par.v;
  p_M = par.p_M;
  k_J = par.k_J;
  del_Me = par.del_Me;
  del_Mt = par.del_Mt;
  fRu = par.fRu;
  w =w;
  d_V=d_V;
  for i=1:length(data.dmetsvl(:,1))
      parss = parsolver;
      %parss.bn=bn;
      parss.Texp = C2K(temparr(i)); % initial calculation to be done at reference temperature then we can retransform
      parss.Tref = T_ref; % reference temperature of AmP values
      TCtemp=tempcorr(parss.Texp, T_ref, [T_A TrefL TrefH TA_l TA_h]);
      pT_Am = TCtemp * p_Am; vT = TCtemp * v; pT_M = TCtemp * p_M; kT_J = TCtemp * k_J; 
      parss.spAm = pT_Am; 
      parss.v = vT; 
      parss.spM = pT_M ; 
      parss.kJ = kT_J;
      t_tmp = [0:0.1:500]';
      parss.terminate=2;
      [to, Xo, ~, phys,TE,YE,IE] = solve_deri(t_tmp, X0, parss);
      durtemp(i)  = TE(3) - TE(2);
      svlstemp(i) = 10. * YE(3,3)/del_Me;
      totltemp(i) = 10. * YE(3,3)/del_Mt;
      masstemp(i) = 1000.*YE(3,3)^3 * d_V * (1 + w * fRu ) / 0.06; % in mg
  end
  % need to introduce them as relative quantities in order account for the
  % different datasets
  prdData.dmetage   = durtemp /durtemp22  ;
  prdData.dmetsvl   = svlstemp/svlstemp22 ;
  prdData.dmettl    = totltemp/totltemp22 ;
  prdData.dmettm    = masstemp/masstemp22 ;



% %% Marini data on AMA tests
fratios = [1, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05]; %feeding treatments wrt to max
X0(1,:) = E_0;     % state 1 is the reserve
X0(2,:) = 0;       % state 3 maturity
X0(3,:) = 1e-5;    % state 4 is the structural length
X0(4,:) = 0;       % state 5 is the cumulative reproduction (NOT USED)
X0(5,:) = 0;       % state 6 is the scaled damage (referenced to water) (NOT USED)
X0(6,:) = 1;       % state 7 is survival probability (NOT USED)
Xkratio = Xkratio;
fMAR = fMAR;
MarAMAsvl_f100=MarAMAsvl_f100;
for fidx=1:numel(fratios)
    parss = parsolver;
    pT_Am = TC_tAMA * p_Am; vT = TC_tAMA * v; pT_M = TC_tAMA * p_M; kT_J = TC_tAMA * k_J; 
    parss.spAm = pT_Am; 
    parss.v = vT; 
    parss.spM = pT_M ; 
    parss.kJ = kT_J;
    epsilon = fMAR - (1 +  Xkratio)^(-1);
    if fratios(fidx) == 1
        newf = fMAR;
    else
        newf = (1 +  Xkratio / fratios(fidx) )^(-1) + epsilon;
    end
    parss.f = [fMAR, newf];


    parss.terminate=0;
    stri = num2str(100*fratios(fidx));
    prnamesvl = ['MarAMAsvl_f',stri];
    if length(MarAMAsvl_f100)>10
        t_tmp = union([0:0.1:37], data.(prnamesvl)(:,1));
    else
        t_tmp = [0:0.1:37];
    end
    [to, Xo, ~, phys,TE,YE,IE] = solve_deriAMA(t_tmp, X0, parss);
    idx = find(ismember(to, data.(prnamesvl)(:,1)));
    svl = 10. .* Xo(idx,3)./del_Me;
    dw = Xo(idx,3).^3 .* d_V .* (1 + w .* phys.e(idx));  % full computation of the weight 
    ww = dw ./ arrayfun(@(x) rescaleweight(x, parsolver.EH42, parsolver.EHj), Xo(idx,2));
    prdsData{fidx} = svl;
    prdwData{fidx} = ww;
end

for fidx=1:numel(fratios)
    stri = num2str(100*fratios(fidx));
    prnamesvl = ['MarAMAsvl_f',stri];
    prnameww = ['MarAMAw_f',stri];
    prdData.(prnamesvl) = prdsData{fidx};
    prdData.(prnameww) = prdwData{fidx};
end

end
