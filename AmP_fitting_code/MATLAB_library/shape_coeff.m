% Author C. Romoli - ibacon GmbH
% 
% change the shape coefficient from that of a tadpole to that of a froglet
function delm = shape_coeff(EH, EH42, EHj, delm1, delm2)
    % delm1 = shape coefficient for tadpole
    % delm2 = shape coefficient for frog
    if (EH<=EH42)
        delm = delm1;
    elseif (EH >= EHj)
        delm = delm2;
    else 
        delm = delm1 + (delm2-delm1)/(EHj - EH42) * (EH - EH42);
    end
end