% Model to analyze the temporal dynamic effects of single stressors 
% Main focus:analysis of effects of chemical substances

function Y = SinglStressor_Model(Cext,parameter,T,N0)

% here example values are applied. You need to adapt the values for your
% pet
clear globals

% parameter values
beta = parameter(1); % effect per time and internal concentration
mu   = parameter(2); % time delay
alfa = parameter(3); % speed of adaptation relative to effect
C   = Cext; % external concentration
param = [beta mu alfa C];

% Anfangswerte
Cint0 = 0;
% NCd0 Inputparameter;
e0 = 0;
a0 = 0;
Y0 = [Cint0;N0;e0;a0];

% Lösung der Differentialgleichungen
[t,Y] = ode45(@(t_t,t_Y) Stressor_ode(t_t,t_Y,param),T,Y0);

%----------------------------------------------------------------------
% Differential equation: dY for dY/dt etc.
function dY = Stressor_ode(t,Y,param)

beta = param(1);
mu = param(2);
alfa = param(3);
C = param(4);

Cint    = Y(1);
N       = Y(2);
e      = Y(3);
a      = Y(4);

%----additional variables needed---%
C0 = 0; % this is the optimum concentration. For hazardous substances this is for example 0,
% for other stressors such as temperature it can have another
% value

% Therefore the zero effect concentraion is set to 0.
CTolext = 5; % e.g. external tolerance concentration in ug/L 
CTolInt = 15; % internal tolerance concentration e.g. ug substance/g mussel weight

%**************************************************************************
% uptake kinetics/ growth dilution
%**************************************************************************

    kup = 0.5; % uptake rate
    kex = 0.03; % excretion rate
    qC = kup/kex; % bioaccumulation factor 
    CTol = CTolInt;% here the internal tolerance concentration is applied

gr_test = 1; % change to 0 if no data are available
if gr_test == 1
    Lc  = 2;% size of the animal 
    Lm  = 15; % maximum lengh of the animal 
    rB  = 0.8; % Bertalanffy growth rate of the animal
    L   = Lm - (Lm-Lc) * exp(-rB *t);
    dL  = rB * (Lm-L); % dL is the change of lenght
    
    dCint   = kup*Lm/L*(C-Cint/qC)- Cint*3/L*dL - kex*Cint;
    
        
else % change if it will be used
   
    dCint   = kupCd*C - kex*Cint;
        
    if Cint > qC %internal concnentration cannot exceed the equilibrium concentration in the body
        Cint = qC;
    end
end


dN      = - eCd*NCd; % eCd should be sth like the percentage of the population
e_equCd   = beta*(Cint-C0)*(Cint>CTol)/(CTol-C0)*(1-aCd);
deCd      = mu*(e_equCd-eCd);
daCd      = (alfa*eCd)*(1-aCd);

dY      = [dCint;dN;deCd;daCd];


