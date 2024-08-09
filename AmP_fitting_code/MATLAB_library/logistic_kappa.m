% Author C. Romoli - ibacon GmbH
%
% function to shift the value of the kap paramter due to stress while still
% keeping it bounded between 0 and 1
function [newkap] = logistic_kappa(kap, stress)
    midpoint = log(kap/(1-kap));
    newkap = 1/(1 + exp(-(stress+midpoint)));
end