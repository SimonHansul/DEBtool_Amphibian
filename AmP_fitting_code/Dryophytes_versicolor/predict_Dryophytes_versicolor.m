function [prdData, info] = predict_Dryophytes_versicolor(par, data, auxData)
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

    % compute temperature correction factors
  TC = tempcorr(temp.am, T_ref, T_A);
  TCRi = tempcorr(temp.Ri, T_ref, T_A);
  TC19 = tempcorr(C2K(19), T_ref, T_A);
  TC22 = tempcorr(C2K(22), T_ref, T_A);
  TC23 = tempcorr(C2K(23), T_ref, T_A);
  TC27 = tempcorr(C2K(27), T_ref, T_A);

  fH= f;
  % hatch
  % initial 
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve  
  E_0 = U_E0 * p_Am;     % J, initial reserve
  
  pars_tpj = [g k l_T v_Hb v_Hj]; % for tadpole
  [t_j, t_b, l_j, l_b, info] = get_tp(pars_tpj, f);

  %% pass the parameter to a new par object
  parsolver.delM1  = 0.1;
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
  parsolver.EH42 = E_H42;     % maturity at start of the metamorph. climax
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

fh=parsolver.f;


%% zero variate data
parsolver.Texp = temp.ah; % initial calculation to be done at reference temperature then we can retransform
  pT_Am = TC * p_Am;
  vT = TC * v;
  pT_M = TC * p_M; 
  kT_J = TC * k_J; 
  parsolver.spAm = pT_Am; 
  parsolver.v = vT; 
  parsolver.spM = pT_M ; 
  parsolver.kJ = kT_J;
  parsolver.f = f;
  parsolver.Tref = T_ref; % reference temperature of AmP values
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
  aT_42 = TE(3);%/ TC;
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


prdData.ah = TE(1);
prdData.ab = TE(2)/TC23*TC; % assume 23 degrees for the age at birth
% assume no temp stress

E_M = pT_Am / vT;
eb = YE(2,1) ./ (E_M * (YE(2,3).^3)) ;
prdData.Wwb = (YE(2,3).^3) .*d_V.* (1+ w * eb ) ./ 0.06;

prdData.Lp = YE(5,3)/del_M;
ep = YE(5,1) ./ (E_M * (YE(5,3).^3)) ;
prdData.Wwp = (YE(5,3).^3) .* (1+ w * ep );


% reproduction
% does not need to be recomputed with the new code, because the
% differences with the solution without the metamorphosis are negligible
% when dealing with the reproduction rate

% life cycle for tadpole
pars_tpj = [g k l_T v_Hb v_Hj]; % for tadpole
[t_j, t_b, l_j, l_b, info] = get_tp(pars_tpj, f);
pars_tp = [g k l_T v_Hb v_Hp]; % for frog
[t_p, t_b, l_p, l_b] = get_tp(pars_tp, f);

% ultimate for toad
l_i = f - l_T;                  % -, scaled ultimate length at f
L_i = L_m * l_i;                % cm, ultimate structural length at f
Lw_i = L_i/ del_M;              % cm, ultimate total length at f
Ww_i = L_i^3 * (1 + f * w);     % g, ultimate wet weight 

pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector at T
RT_i = TCRi * reprod_rate(L_i, f, pars_R);               % #/d, ultimate reproduction rate at T

prdData.Li = Lw_i;
prdData.Ri = RT_i;
prdData.Wwi = Ww_i;


% life span
pars_tm = [g; l_T; h_a/ k_M^2; s_G]; % compose parameter vector at T_ref
t_m = get_tm_s(pars_tm, f, l_b);     % -, scaled mean life span at T_ref
aT_m = t_m/ k_M/ TC;              % d, mean life span at T

prdData.am = aT_m;


% weight at end of metamorphosis
e46 = YE(4,1) ./ (E_M * (YE(4,3).^3)) ;
Ww46 = (YE(4,3).^3) .*d_V.* (1+ w * e46 ) ./ 0.13;
prdData.Wg46 = Ww46;
prdData.W2g46 = Ww46;

ag46 = (TE(4)-TE(2))*TC; % assuming 20C for this
a2g46 = ag46/TC23;
prdData.ag46=ag46;
prdData.a2g46= a2g46;

%% univariate data
% HernRami2018
parsolver.terminate=1;
pT_Am = TC19 * p_Am;
vT = TC19 * v;
pT_M = TC19 * p_M; 
kT_J = TC19 * k_J; 
parsolver.spAm = pT_Am; 
parsolver.v = vT; 
parsolver.spM = pT_M ; 
parsolver.kJ = kT_J;
parsolver.f = fS04;
parsolver.Texp = C2K(19);
parsolver.Tref = T_ref; % reference temperature of AmP values
t_tmp = [0:0.1:1000];
X0(1,:) = E_0;     % state 1 is the reserve
X0(2,:) = 0;      % state 3 maturity
X0(3,:) = 1e-5;   % state 4 is the structural length
X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
[to, Xo, parout, phys,TE,YE] = solve_deri(t_tmp, X0, parsolver);
parsolver.terminate=1;
newt = union([0:0.1:500], data.Ww6w(:,1)+TE(2));
[to, Xo, parout, phys,TE,YE] = solve_deri(newt, X0, parsolver);
idx = find(ismember(to, data.Ww6w(:,1)+TE(2)));
if phys.e==0
    phys.e = 0.*to;
end
dw = Xo(idx,3).^3 .* d_V .* (1 + w .* phys.e(idx));
ww = dw ./ arrayfun(@(x) rescaleweight(x, parsolver.EH42, parsolver.EHj), Xo(idx,2));
prdData.Ww6w = ww;


% beachy1999
parsolver.terminate=1;
pT_Am = TC22 * p_Am;
vT = TC22 * v;
pT_M = TC22 * p_M; 
kT_J = TC22 * k_J; 
parsolver.spAm = pT_Am; 
parsolver.v = vT; 
parsolver.spM = pT_M ; 
parsolver.kJ = kT_J;
parsolver.Texp = C2K(22);
parsolver.Tref = T_ref; % reference temperature of AmP values

% LLL case
parsolver.switch = -1; 
parsolver.f = [fB99L, fB99H];
t_tmp = union([0:0.1:100], data.LLL(:,1)+4); % hatching time directly from paper (it says 4 days)
[to, Xo, parout, phys,TE,YE] = solve_deriBeachy99(t_tmp, X0, parsolver);
idx = find(ismember(to, data.LLL(:,1)+4));
if phys.e==0
    phys.e = 0.*to;
end
dw = Xo(idx,3).^3 .* d_V .* (1 + w .* phys.e(idx));
ww = dw ./ arrayfun(@(x) rescaleweight(x, parsolver.EH42, parsolver.EHj), Xo(idx,2));
prdData.LLL = ww;

parsolver.switch = 1; 
parsolver.f = [fB99L, fB99H];
t_tmp = union([0:0.1:100], data.LHH(:,1)+4); % hatching time directly from paper (it says 4 days)
[to, Xo, parout, phys,TE,YE] = solve_deriBeachy99(t_tmp, X0, parsolver);
idx = find(ismember(to, data.LHH(:,1)+4));
if phys.e==0
    phys.e = 0.*to;
end
dw = Xo(idx,3).^3 .* d_V .* (1 + w .* phys.e(idx));
ww = dw ./ arrayfun(@(x) rescaleweight(x, parsolver.EH42, parsolver.EHj), Xo(idx,2));
prdData.LHH = ww;

parsolver.switch = 0; 
parsolver.f = [fB99L, fB99H];
t_tmp = union([0:0.1:100], data.LLH(:,1)+4); % hatching time directly from paper (it says 4 days)
[to, Xo, parout, phys,TE,YE] = solve_deriBeachy99(t_tmp, X0, parsolver);
idx = find(ismember(to, data.LLH(:,1)+4));
if phys.e==0
    phys.e = 0.*to;
end
dw = Xo(idx,3).^3 .* d_V .* (1 + w .* phys.e(idx));
ww = dw ./ arrayfun(@(x) rescaleweight(x, parsolver.EH42, parsolver.EHj), Xo(idx,2));
prdData.LLH = ww;

% HHH case
parsolver.switch = -1; 
parsolver.f = [fB99H, fB99L];
t_tmp = union([0:0.1:100], data.HHH(:,1)+4); % hatching time directly from paper (it says 4 days)
[to, Xo, parout, phys,TE,YE] = solve_deriBeachy99(t_tmp, X0, parsolver);
idx = find(ismember(to, data.HHH(:,1)+4));
if phys.e==0
    phys.e = 0.*to;
end
dw = Xo(idx,3).^3 .* d_V .* (1 + w .* phys.e(idx));
ww = dw ./ arrayfun(@(x) rescaleweight(x, parsolver.EH42, parsolver.EHj), Xo(idx,2));
prdData.HHH = ww;

parsolver.switch = 1; 
parsolver.f = [fB99H, fB99L];
t_tmp = union([0:0.1:100], data.HLL(:,1)+4); % hatching time directly from paper (it says 4 days)
[to, Xo, parout, phys,TE,YE] = solve_deriBeachy99(t_tmp, X0, parsolver);
idx = find(ismember(to, data.HLL(:,1)+4));
if phys.e==0
    phys.e = 0.*to;
end
dw = Xo(idx,3).^3 .* d_V .* (1 + w .* phys.e(idx));
ww = dw ./ arrayfun(@(x) rescaleweight(x, parsolver.EH42, parsolver.EHj), Xo(idx,2));
prdData.HLL = ww;

parsolver.switch = 0; 
parsolver.f = [fB99H, fB99L];
t_tmp = union([0:0.1:100], data.HHL(:,1)+4); % hatching time directly from paper (it says 4 days)
[to, Xo, parout, phys,TE,YE] = solve_deriBeachy99(t_tmp, X0, parsolver);
idx = find(ismember(to, data.HHL(:,1)+4));
if phys.e==0
    phys.e = 0.*to;
end
dw = Xo(idx,3).^3 .* d_V .* (1 + w .* phys.e(idx));
ww = dw ./ arrayfun(@(x) rescaleweight(x, parsolver.EH42, parsolver.EHj), Xo(idx,2));
prdData.HHL = ww;

% experiment entries High food level (target stage = GS42)
j=1;
dw=[];
for i=1:length(data.temptimeH(:,1))
    parsolver.terminate=1;
    TC=tempcorr(C2K(data.temptimeH(i,1)), T_ref, T_A);
    string = ['timeWwH',num2str(data.temptimeH(i,1))];
    pT_Am = TC * p_Am;
    vT = TC * v;
    pT_M = TC * p_M; 
    kT_J = TC * k_J; 
    parsolver.spAm = pT_Am; 
    parsolver.v = vT; 
    parsolver.spM = pT_M ; 
    parsolver.kJ = kT_J;
    parsolver.f = fH;
    parsolver.Texp = C2K(data.temptimeH(i,1));
    parsolver.Tref = T_ref; % reference temperature of AmP values 
    X0(1,:) = E_0;     % state 1 is the reserve
    X0(2,:) = 0;      % state 3 maturity
    X0(3,:) = 1e-5;   % state 4 is the structural length
    X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
    X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
    X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
    t_tmp = [0:0.1:100];
    [~, ~, ~, ~,TE,~] = solve_deri(t_tmp, X0, parsolver);  % just to get birth timr
    if ismember(data.temptimeH(i,1), [18,20,22,24])
        t_tmp = union(t_tmp, data.(string)(:,1)+TE(2));
    end
    parsolver.terminate=0;
    [to, Xo, ~, ~,TE,YE] = solve_deri(t_tmp, X0, parsolver);
    if ismember(data.temptimeH(i,1), [18,20,22,24])
        To{j} = to;
        Xt{j} = Xo;
        TEb{j} = TE(2);
        j=j+1;
    end
    dw(i) = YE(3,3).^3 .* d_V .* (1 + w .* fH);
    tim(i) = TE(3)-TE(2);
end

Temp=[18,20,22,24];
for i=1:4
    temp = Temp(i);
    string = ['timeWwH',num2str(temp)];
    idx = find(ismember(To{i},data.(string)(:,1) + TEb{i}));
    dwlist = Xt{i}(idx,3).^3 .* d_V .* (1 + w .* fH);
    prdData.(string) = dwlist ./ arrayfun(@(x) rescaleweight(x, parsolver.EH42, parsolver.EHj), Xt{i}(idx,2));
end

for i=1:length(data.temptimeH2(:,1))
    TC=tempcorr(C2K(data.temptimeH2(i,1)), T_ref, T_A);
    pT_Am = TC * p_Am;
    vT = TC * v;
    pT_M = TC * p_M; 
    kT_J = TC * k_J; 
    parsolver.spAm = pT_Am; 
    parsolver.v = vT; 
    parsolver.spM = pT_M ; 
    parsolver.kJ = kT_J;
    parsolver.f = fH;
    parsolver.Texp = C2K(data.temptimeH2(i,1));
    parsolver.Tref = T_ref; % reference temperature of AmP values
    parsolver.terminate=1;
    t_tmp = [0:0.1:1000];
    X0(1,:) = E_0;     % state 1 is the reserve
    X0(2,:) = 0;      % state 3 maturity
    X0(3,:) = 1e-5;   % state 4 is the structural length
    X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
    X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
    X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
    [~, ~, ~, ~,TE,YE] = solve_deri(t_tmp, X0, parsolver);
    dw2(i) = YE(4,3).^3 .* d_V .* (1 + w .* fH);
    tim2(i) = TE(4)-TE(2);
end


prdData.temptimeH = tim';
prdData.tempDwH = dw';
prdData.temptimeH2 = tim2';
prdData.tempDwH2 = dw2';

% experiment entries Low food level (target stage = GS42)
j=1;
for i=1:length(data.temptimeL(:,1))
    parsolver.terminate=1;
    TC=tempcorr(C2K(data.temptimeL(i,1)), T_ref, T_A);
    string = ['timeWwL',num2str(data.temptimeL(i,1))];
    pT_Am = TC * p_Am;
    vT = TC * v;
    pT_M = TC * p_M; 
    kT_J = TC * k_J; 
    parsolver.spAm = pT_Am; 
    parsolver.v = vT; 
    parsolver.spM = pT_M ; 
    parsolver.kJ = kT_J;
    parsolver.f = fL;
    parsolver.Texp = C2K(data.temptimeL(i,1));
    parsolver.Tref = T_ref; % reference temperature of AmP values
    X0(1,:) = E_0;     % state 1 is the reserve
    X0(2,:) = 0;      % state 3 maturity
    X0(3,:) = 1e-5;   % state 4 is the structural length
    X0(4,:) = 0;      % state 5 is the cumulative reproduction (NOT USED)
    X0(5,:) = 0;      % state 6 is the scaled damage (referenced to water) (NOT USED)
    X0(6,:) = 1;      % state 7 is survival probability (NOT USED)
    t_tmp = [0:0.1:100];
    [~, ~, ~, ~,TE,~] = solve_deri(t_tmp, X0, parsolver);
    if ismember(data.temptimeL(i,1), [18,20,22,24])
        t_tmp = union(t_tmp, data.(string)(:,1)+TE(2));
    end
    parsolver.terminate=0;
    [to, Xo, ~, ~,TE,YE] = solve_deri(t_tmp, X0, parsolver);
    if ismember(data.temptimeH(i,1), [18,20,22,24])
        To{j} = to;
        Xt{j} = Xo;
        TEb{j} = TE(2);
        j=j+1;
    end
    dw(i) = YE(4,3).^3 .* d_V .* (1 + w .* fL);
    tim(i) = TE(4)-TE(2);
end
prdData.temptimeL = tim';
prdData.tempDwL = dw';

for i=1:4
    temp = Temp(i);
    string = ['timeWwL',num2str(temp)];
    idx = find(ismember(To{i},data.(string)(:,1)+TEb{i}));
    dwlist = Xt{i}(idx,3).^3 .* d_V .* (1 + w .* fL);
    prdData.(string) = dwlist ./ arrayfun(@(x) rescaleweight(x, parsolver.EH42, parsolver.EHj), Xt{i}(idx,2));
end



end
