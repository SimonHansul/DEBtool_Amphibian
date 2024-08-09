% Author C. Romoli - ibacon GmbH
% 
% to write this code, inspiration from BYOM package by T. Jager
%
% file derivatives.m containing all the derivatives
% this is the core of the program
function dX = derivatives(t, X, par)
    E  = X(1); % state 1 is the reserve
    EH = X(2); % state 2 is the maturity
    L  = X(3); % state 3 is the volumetric length
    % Rc = X(4); % NOT USED state 4 is the cumulative reproduction (not used in this function)
    Dw = X(5); % NOT USED state 5 is the scaled damage (referenced to water)
    S  = X(6); % NOT USED state 6 is survival probability
    
    kapH = 1; % maturation efficiency (set to 100%)
    
    % change kappa and f depending on the maturity stage.
    % This allows change kappa and f somehow arbitrarily during
    % metamorphosis
    kap = par.kap;
    s_nat = - par.bn * (par.Texp - par.Tref)/par.Tref; % calculate stress factor due to temperature
    par.kap2 = logistic_kappa(par.kap, s_nat);         % calculate new kappa value
    if EH >= par.EHb && EH <par.EHj
        % apply new kappa value only after birth. Assumption that embryo is
        % still isolated from stress factor (but still affected by Arrhenious relation)
        kap = par.kap2;
    end
    % function to make sure that we have the right f factor for this
    % developmental stage (f=0 during climax and f=f otherwise)
    [f] = f_function(EH,par.f, par.f2, par.EHb, par.EH42, par.EHj);

    % Calculate energy fluxes
    pA = f * par.spAm * L^2;        % assimilation power
    pS = par.spM * L^3 + par.spT * L^2; % maintenance power (somatic and surface-specific)
    pC = E * (par.EG * par.v * L^2 + pS) / (kap * E + par.EG * L^3); % mobilisation power
    pG = kap * pC - pS;         % growth power
    pJ = par.kJ * EH;           % maturity maintenance power
    pR = (1-kap) * pC - pJ;     % maturation/reproduction power

    % metamorphosis not dealt in this part of the code as we stop the ODE
    % before and transfer the solution to the next code that deals with the
    % metamorphosis

    % Starvation rules may override these fluxes. This section is taken
    % from Tjalling implementation of stdDEB
    if pG < 0        % then we have starvation
        if pC > pS + pJ          % mobilisation is still enough to pay somatic and maturity maintenance
            pG = 0;              % stop growth
            pR = pC - pS - pJ;   % maturation/reproduction gets what's left
        elseif pC >= pS          % mobilisation can only pay for somatic maintenance costs
            pG = 0;              % stop growth
            pR = 0;              % stop maturation/reproduction
            % pJ = pC - pS;        % maturity maintenance gets what's left
        else                     % mobilisation cannot even pay for somatic maintenance
            % pJ = 0;              % stop maturity maintenance
            pR = 0;              % stop maturation/reproduction
            pG = (pC - pS - pJ) / 0.8;%par.yP; % shrinking of structure to pay somatic maintenance and maturity maintenance?
        end
    end
    
    % State variables derivatives
    dE = pA - pC;               % change in reserve
    dL = (1/(3*L^2)) * (pG/par.EG); % change in volumetric length

    if EH < par.EHp  % for embryos and juveniles ...
        dEH = max(0,kapH * pR); % maturation (don't allow to become negative)
        R   = 0;                % no reproduction
    else                        % for adults ...
        dEH = 0;                % no maturation
        R = 0;
        % we do not care about what happens after puberty. The model ends
        % earlier than that
    end
    dRc = R; % change in cumulated reproduction NOT USED
    
    dDw =0; % for now no changes in this as we will add stressors afterwards
    dS =0;
    % update derivatives
    dX = [dE;dEH;dL;dRc;dDw;dS];
end
