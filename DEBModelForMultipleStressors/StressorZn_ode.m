% differential equation for the effect of a single stressor
% dY represents dY/dt etc.
function dY = StressorZn_ode(t,Y,NneyW1Znmod)

% Parameter values - from literature
kexZn = 0.020;          % efflux rate constant after food uptake
                        % + 0.011 efflux rate constant after dissolved
                        % uptake (data by Wang et al 1996, accumulation model)  

kupZn = 1.044;          % uptake rate constant in l g^-1 d^-1, data by Wang et al 1996 
                        % (bioaccumulation model) 
                                               
                        
qCZn = kupZn/kexZn;     % CPbint/CPb) im Gleichgewicht, 
                        % bestimmt kex = kupPb/qCPb, Abbaurate  kex = kup/qC;


CZn0 = 0; % this is the optimum concentration

CZnTolExt = 10; % growth inhibition after 22 days exposure time (Strömgren 1982)

CZnTolInt = 227.1217; % calculated with data above with uptake model

% use the internal tolerance concentraion
CZnTol = CZnTolInt;


% paramters calculated (20.02.2018)
betaZn = 1.309e-05;
muZn = 49.9763;
alfaZn = 9.9867e-07;

% environmental data
CZn = NneyW1Znmod; 


CZnint  = Y(1); 
NZn     = Y(2); 
eZn     = Y(3); 
aZn     = Y(4); 

gr_test = 1; % change to 0 if no data are available
if gr_test == 1
    Lc  = 0.009; % any size, which should be investigated
    % 0.009 = lenght at birth, data from Saraiva et al 2011
    Lm  = 15; % maximum lengh (data from Van der Veer et al 2006)
    rB  = 0.848; % Bertalanffy growth rate (Van der Veer et al 06 parameter estimation)
    L   = Lm - (Lm-Lc) * exp(-rB * t);
    dL  = rB * (Lm-L); % dL is the change of lenght
    dCZnint = kupZn*Lm/L*(CZn-CZnint/qCZn)-CZnint*3/L*dL- kexZn*CZnint;

else
    % uptake kinetics for an individuum
    % g = 0.0848; % Bertalanffy growth rate (Van der Veer et al 06 parameter estimation),
    % Wang and Fischer 1997 used values between 0 and 0.12
    % qCZn = kupZn/kexZn; % concentration at eqilibrium
    % dCZnint   = (kupZn*CZn)/(kexZn+g);%kupCd*(CCd-CCdint/qCCd); %kupCd*(CCd-CCdint/qCCd);
    dCZnint   = kupZn*CZn - kexZn*Znint;
    if CZnint > qCZn %internal concnentration cannot exceed the equilibrium concentration in the body
        CZnint = qCZn;
    end
end


            
%plot(t,L)
%**************************************************************************
% defense in the body (e.g. by lysosomes)
%**************************************************************************



dN      = - eZn*NZn; % eZn should be sth like the percentage the population
                     % size is reduced over time
% e_equZn   = betaZn*(CZnint-CZn0)*(CZnint>CZn0)*(1-aZn);
e_equZn   = betaZn*(CZnint-CZn0)*(CZnint>CZn0)/(CZnTol-CZn0)*(1-aZn);
deZn      = muZn*(e_equZn-eZn);
daZn      = alfaZn*eZn*(1-aZn);

dY      = [dCZnint;dN;deZn;daZn]; 
%----------------------------------------------------------------------

