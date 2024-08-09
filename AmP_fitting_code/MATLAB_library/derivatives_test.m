% Author C. Romoli - ibacon GmbH
% 
% to write this code, inspiration from BYOM package by T. Jager
%
% file derivatives_test.m containing all the derivatives
% This part of the is used to solve from the beginning of the metamorphosis
% climax to the end of the simulation.
function [dX] = derivatives_test(t, X, par)
    E  = X(1); % state 1 is the reserve
    EH = X(2); % state 2 is the maturity
    L  = X(3); % state 3 is the volumetric length
    % Rc = X(4); % NOT USED state 4 is the cumulative reproduction (not used in this function)
    Dw = X(5); % NOT USED state 5 is the scaled damage (referenced to water)
    S  = X(6); % NOT USED state 6 is survival probability
    
    kapH = 1;
    
    % change kappa and f depending on the maturity stage.
    % This allows change kappa and f somehow arbitrarily during
    % metamorphosis
    kap = par.kap;
    s_nat = - par.bn * (par.Texp - par.Tref)/par.Tref;
    par.kap2 = logistic_kappa(par.kap, s_nat);
    if EH >= par.EHb && EH <par.EHj
        kap = par.kap2; % MOD1, simple change of kappa parameter during the whole metamorphosis
    end
    [f] = f_function(EH,par.f, par.f2, par.EHb, par.EH42, par.EHj);
    
    % No starvation rules to be defined here as the model is not so much
    % interested in the stage after the metamorphosis. We assume that
    % afterwards the conditions for the froglet are good and not additional
    % stress should happen

    % State variables derivatives
    % dV1=0; % normal structure
    % dV2=0; % structure to be reprocessed
    % dE2=0;
    dDw = 0; % No damage, NOT USED
    if (EH >= par.EH42 && EH<par.EHj)
        % model with e constant
        % pA is not defined becuse there is no feeding during the met.
        % climax
        pS = par.spM * (L^3); % maintenance power (somatic and surface-specific)
        pC = E * (par.EG * par.v * L^2 + pS) / (kap * E + par.EG * L^3); % mobilisation power
    
        pG = kap * pC - pS;         % growth power
        pJ = par.kJ * EH;           % maturity maintenance power
        pR = (1-kap) * pC - pJ;     % maturation/reproduction power
        
        pL = (E/L^3 * pG/par.EG + pC) / (par.EG * par.kV +E/(L^3)); % pL if e is constant
        % this is a simplification of the expression in the insect article (simplifying everything you get EG there)
        % might become problematic when tox stressors are applied
        dE = pL * par.EG * par.kV - pC;   
        dL = dE * (L/3./E); % this is the equation that would keep e constant
    else 
        % calculation of stress we leave it out for the moment
        pA = f * par.spAm * L^2;        % assimilation power 
        pS = par.spM * L^3 + par.spT * L^2; % maintenance power (somatic and surface-specific)
        pC = E * (par.EG * par.v * L^2 + pS) / (kap * E + par.EG * L^3); % mobilisation power
        pG = kap * pC - pS;         % growth power
        pJ = par.kJ * EH;           % maturity maintenance power
        pR = (1-kap) * pC - pJ;     % maturation/reproduction power
        dL = (1/(3*L^2)) * (pG/par.EG); % change in volumetric length
        dE = pA - pC;               % change in reserve
    end

    % puberty and reproduction..not really important here as we focus on
    % juveniles
    if EH < par.EHp  % for embryos and juveniles ...
        dEH = max(0,kapH * pR); % maturation (don't allow to become negative)
        R   = 0;                % no reproduction
    else                        % for adults ...
        dEH = 0;                % no maturation
        R = 0;
        % we do not care about what happens after puberty. The model ends
        % earlier than that
    end
    dRc = R; % change in cumulated reproduction
    dS =0;
    dX = [dE;dEH;dL;dRc;dDw;dS];
end