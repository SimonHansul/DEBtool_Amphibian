% Author C. Romoli - ibacon GmbH
% 
% Function to rescale wet weight accounting for a changing ratio dry/wet
% weight along the metamorphosis climax
% Arguments:
% - EH : maturity level
% - EH42 : maturity at start of the metamorphosis
% - EHj : maturity at end of the metamorphosis
function DW_ratio = rescaleweight(EH,EH42,EHj)
    if EH<EH42
        DW_ratio = 0.06;
    elseif EH >= EH42 && EH<EHj
        DW_ratio = (0.13 - 0.06)/(EHj-EH42) * (EH-EH42)+0.06;
    else 
        DW_ratio = 0.13;
    end
end