% differential equation for the effect of a single stressor
% dY represents dY/dt etc.
function dY = StressorCu_ode(t,Y,NneyW1Cumod)


    kexCu = 0.0990; % elimination rate under equilibrium conditions
    % calculated with the formula log(2)/biological half life of Cu from
    % data by Adema 1981 in Haren et al 1990. However, model performance of
    % Haren not satisfying
    % other option: 1.0560; from JI Lorenzo, E Aierbe, VK Mubiana, R Blust, R Beiras
% Molluscan Shellfish Safety, 533-544, the value was corrected for the
% correct unit. Unfortunately this value
% overestimates the elimination rate as it can be seen in the same
% publication in a comparison of their model with observational data so
% both options are not perfect
    kupCu =  0.2921; % uptake rate under equilibrium conditions
    % calculated with the formula bioconcentration factor (BC) = uptake rate/elimination rate from
    % --> uptake rate = BC*kex, data by Adema 1981 in Haren et al 1990.
    % other option: 0.6960;  % from JI Lorenzo, E Aierbe, VK Mubiana, R Blust, R Beiras
qCCu = 2.95; % =kupCu/kexCu;% is  2.95 per 1 g adw (data found in Haren et al 1990, 
% original data (in dutch) by Adema 1981 (Accumulatie  en  eliminatie  van 
% enkele metalen door  de mossel Mytilus  edulis  laboratorium onderzoek, SD 0.3 
% kex = kup/qC; equilibrium, kex is excretion 

CCu0 = 0; % this is the optimum concentration. For Cd this is for example 0,
% for other potential stressors such as temperature it can have another
% value

% Therefore the 0 effect concentraion is set to 0.
CCuTolext = 3; % in ug/L - tolerance value based on Strömgren 1982
% (growth inhibition after 12 days exposure)

% modelled internal concentration based on data
% by Strömgren and modelled with DEB uptake model (same structure as below,
% size and growth dilution considered
CCuTolInt = 2.37354694729913;


    CCuTol = CCuTolInt;%*qCCd; % here the bioaccumulation factor is used to translate the 
    % external concentration in an internal concentration.


% calculated parameters (4.6.2018)
betaCu = 0.051964121197595;            % mortality rate per time and internal conc
muCu = 5.589424382286022e-04;              % 1/latency time
alfaCu = 1.733928496063473e-05;        % acclimation rate, relative to observation


% environmental data
CCu = NneyW1Cumod; 


CCuint  = Y(1); 
NCu     = Y(2); 
eCu     = Y(3); 
aCu     = Y(4); 


% model
gr_test = 1; % change to 0 if no data are available
if gr_test == 1
    Lc  = 0.009; % any size, which should be investigated
    % 0.009 = lenght at birth, data from Saraiva et al 2011
    %L0 (lenght at start of development) from Daphnia here changed to Lb (lenght at birth 
    % (means start of feeding)

    Lm  = 15; % maximum lengh (data from Van der Veer et al 2006)
    rB  = 0.848; % Bertalanffy growth rate (Van der Veer et al 06 parameter estimation)
    L   = Lm - (Lm-Lc) * exp(-rB * t); 
    dL  = rB * (Lm-L); % dL is the change of lenght
    
    % This is the size corrected uptake rate minus the growth dilution and
    % minus the elimination rate
    dCCuint = kupCu*Lm/L*(CCu-CCuint/qCCu)-CCuint*3/L*dL- kexCu*CCuint;
    
else
% uptake kinetics for an individuum

    g = 0.0848; % Bertalanffy growth rate (Van der Veer et al 06 parameter estimation),
    % Wang and Fischer 1997 used values between 0 and 0.12
    
    dCCuint   = kupCu*CCu - kexCu*Cuint;
        
    if CCuint > qCCu %internal concnentration cannot exceed the equilibrium concentration in the body
        CCuint = qCCu;
    end
end

%**************************************************************************
% changes 
%**************************************************************************

dN      = - eCu*NCu;
e_equCu   = betaCu*(CCuint-CCu0)*(CCuint>CCu0)/(CCuTol-CCu0)*(1-aCu);
deCu      = muCu*(e_equCu-eCu);
daCu      = alfaCu*eCu*(1-aCu);

dY      = [dCCuint;dN;deCu;daCu]; 
%----------------------------------------------------------------------

