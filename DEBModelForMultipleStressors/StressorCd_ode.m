% differential equation for the effect of a single stressor
% dY represents dY/dt etc.
function dY = StressorCd_ode(t,Y,StationCdmod)

% parameters from literature data
    kupCd = 0.578; % from Wang and Fisher 97 (value for 2.5 cm size 
    % mussels = 0.578 in g/ day)
    kexCd = 0.0324; % per day % (value for 2.5 cm size mussels 0.0324) from Wang and Fisher 1997
    % uptake kinetics for an individuum
    qCCd = kupCd/kexCd; % Bioaccumulation factor 

    CCd0 = 0; % this is the optimum concentration. For Cd this is for example 0,
    % for other potential stressors such as temperature it can have another
    % value
    
CCdTolInt = 16.9134486961835;
CCdTol = CCdTolInt;
 % ug/g mussel weight % calcluated from the external tolerance
% concentration with the script SchadstoffCalcExtIntTol and based on data
% by Wang and Fischer 1997 for mussels of the size 2.25cm

% calculated parameters (19.02.2018)
betaCd = 0.001634242590277;% 0.1589;         % mortality rate per time and internal conc
muCd = 0.005054977743254;                    % 1/latency time
alfaCd = 0.999932705204765;                  % acclimation rate, relative to mortality rate


% environmental data
w = 1;
if w == 1   % if monitoring data about the metal concentration are available 
CCd = StationCdmod; % use these data
else      % if such data are not available use the sediment data and predict 
          % the environmental concentration with the help of the water 
          % column partition coefficient kd
CCd = 0;
disp('values for concentrations of Cd missing')
end
          

CCdint  = Y(1);
NCd     = Y(2);
eCd     = Y(3);
aCd     = Y(4);
gr_test = 1; % change to 0 if no data are available
if gr_test == 1
    Lc  = 0.009; % any size, which should be investigated
    % 0.009 = lenght at birth, data from Saraiva et al 2011
    %L0 (lenght at start of development) from Daphnia here changed to Lb (lenght at birth
    % (means start of feeding)
    % Lp  = 1.2; % lenght at puberty (data from Saraiva et al 2011)
    Lm  = 15; % maximum lengh (data from Van der Veer et al 2006)
    rB  = 0.848; % Bertalanffy growth rate (Van der Veer et al 06 parameter estimation)
    L   = Lm - (Lm-Lc) * exp(-rB *t);
    dL  = rB * (Lm-L); % dL is the change of lenght

    
    dCCdint   = kupCd* Lm/L* (CCd-CCdint)  - CCdint  * 3/L * dL- kexCd*CCdint;

    % equation part "CCdint*3/L *dl" for growth dilution from
    % implementation from DEBtoxModel original (unter Programme DeBtox - engine - derivates) 
 
    % da dies kein Experiment ist, in dem die Konzentration eines kleinen
    % Volumens sich rasch verringert, sondern hier das Verhalten in einem
    % Meer widergespiegelt werden soll wurde die ursprüngliche Gleichung
    % von Kooijmann mit dieser unter der Annahme, dass sich die 
    % Konzentration in Wasser durch die Aufnahme des Organismus nicht 
    % wesentlich verringert ersetzt.
    % Allerdings wird hier berücksichtigt, dass eine gewisse Menge durch
    % die Ausscheidung etc. wieder verloren geht (letzer Term).
       
else
    kexCd = 0.0324; % (value for 2.5 cm size mussels 0.0324) from Wang and Fisher 1997
    % uptake kinetics for an individuum
    kupCd = 0.578;  % from Wang and Fisher 97 (value for 2.5 cm size mussels = 0.578)
    g = 0.0848; % Bertalanffy growth rate (Van der Veer et al 06 parameter estimation),
    % Wang and Fischer 1997 used values between 0 and 0.12
    qCCd = kupCd/kexCd; % concentration at eqilibrium % Bioaccumulation factor
    
    dCCdint   = kupCd*CCd - kexCd*CCdint;
        
    % old by Ebenhöh: dCCdint = kup*(C-Cint/qC);
    if CCdint > qCCd %internal concnentration cannot exceed the equilibrium concentration in the body
        CCdint = qCCd;
    end
end
%plot(t,L)

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
