
function Y = Stressor_modelCd_(CCdext,parameter,T,NCd0)

clear globals
% Antrieb:  CCdext = externe Schadstoffkonzentration, unten einfach CCd
% parameter = Parameterliste s.u.
% T = Liste von Zeitpunkten, NCd0 = Anfangspopulation
% Y = [CCdint,N,m,a], Variablenwerte an den Zeitpunkten T (Ergebnis)

% Variable: CCdint = interne Schadstffkonzentration, Anfangswert 0
%           N = Populationsgöße, Anfangswert 100
%           m = Mortalität durch Schadstoff (Funktion von CCdint und a)
%           a = Anpassung an den Schadstoff, Anfangswert 0, maximal 1

% Das Modell besteht aus den Teilmodellen
% (1) Schadstoffdynamik:    d/dt(CCdint) = kupCd*CCd - kex*CCdint
% (2) Populationsdynamik:   d/dt(N) = - m*N   Regeneration vernachlässigt
% (3) Physiologische Dynamik:
%             Mortalität:   mort = betaCd*(CCdint-CCd0)*(CCdint>CCd0)*(1-a)
%                           d/dt(m) = muCd*(mort-m)
%             Adaption:     d/dt(a) = alfaCd*m*(1-a)
%             "Zeigervariable": Stress = (CCdint-CCd0)*(CCdint>CCd0)/CCdT
%             CCdT = Toleranzkonzentration, kein echter Parameter
%                           mort = (betaCd*CCdT)*Stress*(1-a)
%                           d/dt(a) = (alfaCd*betaCd*CCdT)*Stress*(1-a)

% Erklärung:
%     Der Stress ist eine Funktion der int. Schadstoffkonz. CCdint
%     Mortalität hängt vom Stress ab, reduziert durch Adaption
%     Verzögerung der Wirkung bei muCd<unendlich, sonst Sofortwirkung m=mort
%     Die Anpassung a entwickelt sich bei Mortalität durch Schadstoff

% Parameter:    kupCd, qCCd, CCd0, betaCd, muCd, alfaCd
% kupCd   = Aufnahmerate
% qCCd    = (CCdint/CCd) im Gleichgewicht, bestimmt kex = kupCd/qCCd, Abbaurate
% CCd0    = Schwellenwert (von CCdint) für Giftwirkung
% CCdT    = Toleranzkonzentration, in betaCd enthalten
% betaCd  = Mortalitätsrate pro CCdT, Dimension 1/(Zeit*Konz)
% muCd    = 1/Latenzzeit, kann unendlich sein, dann m=mort
% alfaCd  = Adaptionsrate relativ zur Mortalität (dimensionslos)

% Parameterwerte

% by exponential function % unit:liter/g/day
% qCCd   = parameter(1);% if no information on the efflux rate is available: kex = kupCd/qCCd;
betaCd = parameter(1);               % pro Zeit und int. Konz.
muCd   = parameter(2);
alfaCd = parameter(3);               % Adaptionsgeschw. rel. zur Mortalität
CCd    = CCdext;
param = [betaCd muCd alfaCd CCd];

% Anfangswerte
CCdint0 = 0;
% NCd0 Inputparameter;
eCd0 = 0;
aCd0 = 0;
Y0 = [CCdint0;NCd0;eCd0;aCd0];

% Lösung der Differentialgleichungen
[t,Y] = ode45(@(t_t,t_Y) StressorCd_ode(t_t,t_Y,param),T,Y0);
%[t,Y] = ode45(@SchadstoffCd_ode,T,Y0);

%----------------------------------------------------------------------
% Differentialgleichung: dY steht für dY/dt usw.
function dY = StressorCd_ode(t,Y,param)
%global CCd betaCd muCd alfaCd L

betaCd = param(1);
muCd = param(2);
alfaCd = param(3);
CCd = param(4);

CCdint    = Y(1);
NCd       = Y(2);
eCd       = Y(3);
aCd       = Y(4);

%----additional variables needed---%
CCd0 = 0; % this is the optimum concentration. For Cd this is for example 0,
% for other potential stressors such as temperature it can have another
% value

% Therefore the 0 effect concentraion is set to 0.
CCdTolext = 5; % in ug/L - tolerance value based on Strömgren 1982
% (growth inhibition after 9 days exposure)
% pCCdTolInt = 18.1226;
CCdTolInt = 16.9134486961835;
 % ug/g mussel weight % calcluated from the external tolerance
% concentration with the script SchadstoffCalcExtIntTol and based on data
% by Wang and Fischer 1997 for mussels of the size 2.25cm

%**************************************************************************
% degradation
%**************************************************************************
% kd = 0.01; % degradation rate constant (use 0 to bypass); default from GUT model
%
% if kd > 0,
%    CCd = CCd * exp(-kd*t); % first order degradation
% end

%**************************************************************************
% uptake kinetics/ growth dilution
%**************************************************************************

    kupCd = 0.578; % from Wang and Fisher 97 (value for 2.5 cm size 
    % mussels = 0.578 in g/ day)
    kexCd = 0.0324; % per day % (value for 2.5 cm size mussels 0.0324) from Wang and Fisher 1997
    % uptake kinetics for an individuum
    qCCd = kupCd/kexCd; % Bioaccumulation factor (Umwandlungfaktor für 
    % interne und externe Konzentration (Conc im Körper = qCCd* conc im
    % Wasser)
    CCdTol = CCdTolInt;%*qCCd; % here the bioaccumulation factor is used to translate the 
    % external concentration in an internal concentration.

gr_test = 1; % change to 0 if no data are available
if gr_test == 1
    Lc  = 2.25;% size of mussels of Publication Súnila et al any size, which should be investigated
    % 0.009 = lenght at birth, data from Saraiva et al 2011
    %L0 (lenght at start of development) from Daphnia here changed to Lb (lenght at birth
    % (means start of feeding)
    % Lp  = 1.2; % lenght at puberty (data from Saraiva et al 2011)
    Lm  = 15; % maximum lengh (data from Van der Veer et al 2006)
    rB  = 0.848; % Bertalanffy growth rate (Van der Veer et al 06 parameter estimation)
    L   = Lm - (Lm-Lc) * exp(-rB *t);
    dL  = rB * (Lm-L); % dL is the change of lenght
    %kupCd  = 2.6056*L^-1.42; % ; %update Dez 2017 equation derived
    % from Wang and Fisher 97 (value for 2.5 cm size mussels = 0.578)
    %kexCd = -0.0012*L^2+0.0027*L; %metal efflux rate constant d^-1 (value 
    % for 2.5 cm size mussels 0.0324); % Wang and Fischer 1997

    kexCd = 0.0324; % (value for 2.5 cm size mussels 0.0324) from Wang and Fisher 1997
    % uptake kinetics for an individuum
    kupCd = 0.578; % from Wang and Fisher 97 (value for 2.5 cm size 
    % mussels = 0.578 in g/ day)
    
    qCCd = kupCd/kexCd; % Bioaccumulation factor (Umwandlungfaktor für 
    % interne und externe Konzentration (Conc im Körper = qCCd* conc im
    % Wasser)

    dCCdint   = kupCd*Lm/L*(CCd-CCdint/qCCd)- CCdint*3/L*dL - kexCd*CCdint;
    
    % da dies kein Experiment ist, in dem die Konzentration eines kleinen
    % Volumens sich rasch verringert, sondern hier das Verhalten in einem
    % Meer widergespiegelt werden soll wurde die ursprüngliche Gleichung
    % von Kooijmann (?) mit dieser unter der Annahme, dass sich die 
    % Konzentration in Wasser durch die Aufnahme des Organismus nicht 
    % wesentlich verringert ersetzt.
    % Allerdings wird hier berücksichtigt, dass eine gewisse Menge durch
    % die Ausscheidung etc. wieder verloren geht (letzer Term).
    
else % 
    kexCd = 0.0324; % (value for 2.5 cm size mussels 0.0324) from Wang and Fisher 1997
    % uptake kinetics for an individuum
    kupCd = 0.578;  % from Wang and Fisher 97 (value for 2.5 cm size mussels = 0.578)
    g = 0.0848; % Bertalanffy growth rate (Van der Veer et al 06 parameter estimation),
    % Wang and Fischer 1997 used values between 0 and 0.12
    qCCd = kupCd/kexCd; % concentration at eqilibrium % Bioaccumulation factor (Umwandlungfaktor für 
    % interne und externe Konzentration (Conc im Körper = qCCd* conc im
    % Wasser)

    dCCdint   = kupCd*CCd - kexCd*CCdint;

    if CCdint > qCCd %internal concnentration cannot exceed the equilibrium concentration in the body
        CCdint = qCCd;
    end
end


dN      = - eCd*NCd; % eCd should be sth like the percentage of the population
e_equCd   = betaCd*(CCdint-CCd0)*(CCdint>CCdTol)/(CCdTol-CCd0)*(1-aCd);
deCd      = muCd*(e_equCd-eCd);
daCd      = (alfaCd*eCd)*(1-aCd);

dY      = [dCCdint;dN;deCd;daCd];
%----------------------------------------------------------------------

