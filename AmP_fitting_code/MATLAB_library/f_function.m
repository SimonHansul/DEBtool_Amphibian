% Author C. Romoli - ibacon GmbH
% 
% codes the lack of feeding furing the climax
function [newf] = f_function(EH, f1,f2,EH1,EH2,EHj)
    if EH < EH1
        newf = 0; % before birth
    elseif (EH < EH2 && EH >= EH1)
        newf = f1; % between birth and start of climax
    elseif EH <EHj && EH >= EH2
        newf = f2; % during climax
    else
        newf = f1; % after climax
    end
end
