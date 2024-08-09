% Author C. Romoli - ibacon GmbH
% 
% change kappa and depending on the maturation stage
% This function was inspired by the work on Crinia georgiana from Mueller
% et al., 2012.
% arguments:
% - EH: maturity
% - kappa1: value of kappa and f before maturity EH1
% - kappa2: value of kappa and f between maturity EH1 and EHj
% - EH1 and EHj: intermediate maturity thresholds (normally they could be the
% start and stop of hte climax)
% - EHj : maturity threshold of end of the metamorphosis
function [newkappa] = kappa_function(EH,kappa1,kappa2,EH1,EHj)
    if (EH >= EH1) && (EH < EHj)
        newkappa = kappa2;
    else
        newkappa = kappa1; % returns to normal level immediately
    end
end

