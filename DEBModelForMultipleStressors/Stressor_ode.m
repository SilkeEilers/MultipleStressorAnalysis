% differential equation for the effect of a single stressor
% dY represents dY/dt etc.
function dY = Stressor_ode(t,Y,EnvConc)

% parameters from literature data
    kup = 0.5; % uptake rate
    kex = 0.03; % excretion rate

    CCd= 0; % this is the optimum concentration. For hazardous substances this is for example 0,
    % for other potential stressors such as temperature it can have another
    % value
    
CTolInt = 16; % internal tolerance value
CTol = CTolInt; % external concentration, tolerance value

% calculated parameters % example
beta = 0.0001;   % mortality rate per time and internal conc
mu = 0.005;      % 1/latency time
alfa = 0.9;    % acclimation rate

% environmental data
w = 1;
if w == 1   % if monitoring data about the metal concentration are available 
C = EnvConc; % ennvironmental concentration
else      % if such data are not available use the sediment data and predict 
          % the environmental concentration with the help of the water 
          % column partition coefficient kd
C = 0;
disp('values for concentrations of stressor missing')
end
          

Cint  = Y(1);
N     = Y(2);
e     = Y(3);
a     = Y(4);
gr_test = 1; % change to 0 if no data are available
if gr_test == 1
    Lc  = 0.009; % any size, which should be investigated
    Lm  = 15; % maximum lengh of the animal
    rB  = 0.8; % Bertalanffy growth rate 
    L   = Lm - (Lm-Lc) * exp(-rB *t);
    dL  = rB * (Lm-L); % dL is the change of lenght
    
    dCint   = kup* Lm/L* (C-Cint)  - Cint  * 3/L * dL- kex*Cint;
       
else

    qCCd = kupCd/kexCd; % bioaccumulation factor 

    dCCdint   = kupCd*CCd - kexCd*CCdint;

    if CCdint > qCCd %internal concnentration cannot exceed the equilibrium concentration in the body
        CCdint = qCCd;
    end
end

dN      = - eCd*NCd; % eCd should be sth like the percentage the population
e_equCd   = betaCd*(CCdint-CCd0)*(CCdint>CCd0)/(CCdTol-CCd0)*(1-aCd);
deCd      = muCd*(e_equCd-eCd);
daCd      = (alfaCd*eCd)*(1-aCd);

dY      = [dCCdint;dN;deCd;daCd];
%------------------------


% The uptake of substances in the organism is calculated so that
% the model taking ino account the bioconcentrationfactor and kinetics
% can be compared to the environmental data; the uptake of the toxicant is
% dependent on exposure time and lenght, other variables are dependent on
% the type of substance and the model organism but are not dynamic


%----------------------------------------------------------------------
