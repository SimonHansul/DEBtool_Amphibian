% function to address the feeding regimes for teh AMA tests performed by
% J.Marini in the paper DOI: 10.1002/etc.5596
function [f_ama] = f_amatest(t, f, fmax)
    if t>15
        f_ama = f;
    else
        f_ama = fmax;
    end
end