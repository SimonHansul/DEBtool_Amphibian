% Author C. Romoli - ibacon GmbH
% 
% to write this code, inspiration from BYOM package by T. Jager
function [tout, Xout, par, phys, TE, YE, IE] = solve_deriAMA(t, X0, par)
    % function to call the ODE solver for the amphibians.
    % input parameters are:
    % - t  : time array
    % - X0 : initial conditions
    % - par: DEB paramters
    %
    % Returned quantities:
    % tout  : time array of the solution 
    % Xout  : array of the state variables of the solution 
    % par   : updated struct with the DEB paramters
    % phys  : structure with some derived quantities (e.g., e the scaled reserve)
    % TE    : array with the events (hatch, birth, beg. of climax, end of metamorphosis)
    % YE    : array with the state variable values corresponding to each event (hatch, birth, beg. of climax, end of metamorphosis)
    % IE    : array of event identifier


    % set up the ODE solver settings and solve the first part of the
    % development (until climax is reached)
    E_M = par.spAm / par.v;  %max reserve
    
    options = odeset(Events=@reach_climax, RelTol=1e-5, AbsTol=1e-5);
    [tout_tmp1,Xout_tmp1,TE,YE,IE] = ode45(@derivativesAMA,t,X0,options,par);

    Xout_tmp = Xout_tmp1;

    % compute wet weight
    E_M = par.spAm / par.v;  %max reserve
    wE = 23.9;
    muE = 550e3;
    dE = par.dV;
    omegaV = E_M * wE / (dE * muE);
    e = Xout_tmp(:,1) ./ (E_M * (Xout_tmp(:,3).^3)) ;
    wwt = (Xout_tmp(:,3).^3) .* (1+omegaV*e);
    phys.e = e;
    % convert structural length in physical length and save teh results into
    % a new structure
    L = Xout_tmp(:,3);
    EH = Xout_tmp(:,2);

     if YE(end,2)<par.EH42 
        tout = tout_tmp1;
        Xout = Xout_tmp1;
        TE(3:4)=[100 100];
        YE(3:4,:) = repelem([10000 10000 10000 10000 10000 10000],2,1);
        IE(3)=3;
        IE(4)=4;
        return;
     end
    % stop the calculation when we reach the start of the metamorphosis
    % climax. This is needed for cases in which the data tell quite clearly
    % that we cannot reach the climax too early and continnue the solution
    % post-climax
     if par.terminate==2
         % avoid bad solution in AmP and stop calculation when entering
         % climax.
         % If statements below needed to account for possibly different
         % shapes of the time arrays
         % for the eventual solution after the climax unreasonable values
         % are given
         if (size(t,1)==1)
             tout = [tout_tmp1(1:end-1); t(t>tout_tmp1(end))'];
         else
             tout = [tout_tmp1(1:end-1); t(t>tout_tmp1(end))];
         end
         Xout = [Xout_tmp1(1:end-1,:) ; repelem([10000 10000 10000 10000 10000 10000],length(t(t>tout_tmp1(end))),1)];
         %phys=0;
         tha=TE(IE==1);
         Yha=YE(IE==1,:);
         if length(tha)<1
             tha=0;
             Yha=X0(1:6)';
         end
         tbi=TE(IE==2);
         t42=TE(IE==3);
         Ybi = YE(IE==2,:);
         Y42 = YE(IE==3,:);
         TE = [tha; tbi; t42];
         YE = [Yha; Ybi; Y42];
         return;
     end

    %% Second part of the solution to account for the metamorphosis
    % Here we start from the final point of reaching the metamorphosis climax
    % fixes different array sizes
    if size(t,1)==1
        t = [0 (t(t>tout_tmp1(end))-tout_tmp1(end))];
    else
        t = [0; (t(t>tout_tmp1(end))-tout_tmp1(end))]';
    end 
    X0 = Xout_tmp1(end,:);
    options = odeset(Events=@end_metamorph, RelTol=1e-5, AbsTol=1e-5);
    [tout_tmp2,Xout_tmp2,TE2,YE2,IE2] = ode45(@derivatives_testAMA,t,X0,options,par);
    if length(t)==2
        % avoid problems when only 2 time points are given
        tout_tmp2 = [tout_tmp2(1) tout_tmp2(end)];
        Xout_tmp2 = [Xout_tmp2(1,:); Xout_tmp2(end,:)];
    end
    
    tout = vertcat(tout_tmp1(1:end-1), tout_tmp2(2:end)+tout_tmp1(end));
    Xout = vertcat(Xout_tmp1(1:end-1,:), Xout_tmp2(2:end,1:6));

    e = Xout(:,1) ./ (E_M * (Xout(:,3).^3));
    
    % new attributes of the structure
    phys.e = e;


if ~isempty(TE2)
    TE(4) = TE2(1)+tout_tmp1(end);
    YE(4,:) = YE2(1,1:6);
    IE(4) = 4;
    if length(TE2)>1
        TE(5) = TE2(end)+tout_tmp1(end);
        YE(5,:) = YE2(2,1:6);
        IE(5) = 5;
    else
        TE(5) = 2000;
        YE(5,:) = [20000 20000 20000 20000 20000 20000];
        IE(5) = -1;
    end
else
    TE(4)=1000;
    YE(4,:) = [10000 10000 10000 10000 10000 10000];
    IE(4)=0;
    TE(5) = 2000;
    YE(5,:) = [20000 20000 20000 20000 20000 20000];
    IE(5) = -1;
end
tha=TE(IE==1);
Yha=YE(IE==1,:);
if length(tha)<1
    tha=0;
    Yha=X0(1:6);
end
tbi=TE(IE==2);
t42=TE(IE==3);
t46=TE(4);
tp =TE(5);

Ybi = YE(IE==2,:);
Y42 = YE(IE==3,:);
Y46 = YE(4,:);
Yp  = YE(5,:);
TE = [tha; tbi; t42; t46; tp];
YE = [Yha; Ybi; Y42; Y46; Yp];
end


% event functions
function [value,isterminal,direction] = end_metamorph(t,X,par)
    EHj  = par.EHj;
    EHp  = par.EHp;
    nevents = 2;
    value = zeros(nevents,1);    
    value(1) = X(2) - EHj;
    value(2) = X(2) - EHp;
    direction = zeros(nevents,1);
    isterminal = [par.terminate, 1];%zeros(nevents,1);
end

function [value,isterminal,direction] = reach_climax(t,X,par)
    EHh  = par.EHh; % hatch
    EHb  = par.EHb;
    EH42 = par.EH42;
    nevents = 3;
    value = zeros(nevents,1);    
    value(1) = X(2) - EHh;
    value(2) = X(2) - EHb;
    value(3) = X(2) - EH42;
    direction = zeros(nevents,1);
    isterminal = zeros(nevents,1);
    isterminal(3,1) = 1; % terminate the ODE once beginning of climax is reached
end
