% differential equation for the effect of a single stressor
% dY represents dY/dt etc.
function dY = StressorPb_ode(t,Y,NneyW1Pbmod,TactCelsius)

kupPb = 0.132;          % 0.132; paper by Schulz-Baldes 1974        
                        % other possibility: 0.515; data by Mubiana et al 2007 
                        % mean of uptake rate at 6°C and 16°C 
                                               
kexPb = 0.013;          %  0.013; elimination/ excretion/ loss
                        % excetion rate constant of paper by Schulz-Baldes 1974 
                        % other option: 0.0031: kex calculated with uptake rate by Schulz-Baldes and
                        % the bioconcentration factor from Han Zhao-Xiang
                        % et al 2013
                        % other possibility: 0.059; data by Mubiana et al 2007 (pooled data, no 
                        % temperature effect
                        % mean of uptake rate at 6°C and 16°C                  
                        % other possibility 2: 0.013; elimination/ excretion/ loss
                        % excetion rate constant paper by Schulz-Baldes 1974 
                        
qCPb = kupPb/kexPb;     % = ca. 10,1538 (CPbint/CPb) im Gleichgewicht, 
                        % bestimmt kex = kupPb/qCPb, Abbaurate  kex = kup/qC;

 
CPbTolext = 0.5; % 30; % in ug/L - reaction after 7 days exposition, before that 
% no statistically significant mortality observed (Zhao-Xiang Han et al.
% 2013)

CPb0 = 0; % this is the optimum concentration

CPbTolInt = 8.157467532467532; % value taken directly from graph in publicatiom
% in ug/g - reaction after 7 days exposition (mortality,
% reduced zymosan phagozytose activity (particles/ mg)(Zhao-Xiang Han et al.
% 2013, Fig. 2)%


CPbTol = CPbTolInt;


% calculated parameters (04.06.2018)
betaPb = 1.231747852031546e-04;
muPb = 8.555054814402578e-04;
alfaPb = 0.999909505701230;

% environmental data
CPb= NneyW1Pbmod; 


CPbint  = Y(1); 
NPb     = Y(2); 
ePb     = Y(3); 
aPb     = Y(4); 

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
    
    % without Temperature interaction
    %dCPbint = kupPb*Lm/L*(CPb-CPbint/qCPb)-CPbint*3/L*dL- kexPb*CPbint;
    
    % including Temperature interaction (here the interaction during uptake
    % process) (new)
    dCPbint = kupPb*Lm/L*(CPb-CPbint/qCPb)+ TactCelsius*0.0018...
        - CPbint*3/L*dL- kexPb*CPbint;
    
else
% uptake kinetics for an individuum
    % dCPbint   = kupPb*(CPb-CPbint/qCPb);
    dCPbint   = kupPb*CPb - kexPb*Pbint;
        
    if CPbint > qCPb %internal concnentration cannot exceed the equilibrium concentration in the body
        CPbint = qCPb;
    end
end


dN      = - ePb*NPb; % ePb should be sth like the percentage the population
                     % size is reduced over time
e_equPb   = betaPb*(CPbint-CPb0)*(CPbint>CPb0)/(CPbTol-CPb0)*(1-aPb);
dePb      = muPb*(e_equPb-ePb);
daPb      = alfaPb*ePb*(1-aPb);

dY      = [dCPbint;dN;dePb;daPb]; 
%---------------------------------------------
%----------------------------------------------------------------------