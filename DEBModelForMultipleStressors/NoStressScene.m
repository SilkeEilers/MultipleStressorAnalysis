% (system of differential equations for basic physiological processes and
% processes over the life span of an individual)
% outputs

function [MV_glob,ME_glob,MR_glob,MH_glob,L] = NoStressScene(tspan, tstart, tend, Tact, Stationphyto)
format compact

global globCRControl
global globCRControlTime
global globJxiI

globJxiI = [];
globCRControl = [];
globCRControlTime =[];

%--------------------------------------------------------
% Temperature correction for the control scenario

% if the actual temperature is greater than the tolerance value, then the
% tolerance temperature is used in this control szenario as above this
% value negative effects would occur

% delete all values above 25°C for the control scenario
Tact(Tact>298)=298;
TactCelsius = double(Tact)-273;
% startDate = datenum('21-02-2005');
% endDate = datenum('21-11-2005');
% xData = linspace(startDate, endDate, tspan);
% figure
% plot(xData,Tact)

figure(1)
plot(tspan,TactCelsius)
title('temperature in control scenario')
xlabel('time in days')
ylabel('temperature in °C')
%--------------------------------------------------------

%==========================================================================
% The organism level - DEB model (Saraiva el al 2012)
%==========================================================================

%==========================================================================
% initial conditions
%==========================================================================
%
% initial conditions for the processes, which should be calculated (at
% birth)/ state variables

disp('DEB model')

MV0 = 3.3e-09; % structural mass at birth in mol C data van der Veer et al 2006;%8.86388638664122e-11; % JVG0;

ME0 = 1.48e-10; % initial reserve mass at optimal food conditions in mol C^E

MH0 = 4.2898e-11; % cumulative maturity at birth, calculated with data 
                  % from Saraiva et al. 2012  Ehb/uE

dv = 0.2; % g dw/ cm^3; Rosland et al 09 in Saraiva et al 2011 dv = dE, 
          % bivalve structure and reserve specific density

wv = 25.22; % f dw/mol, bivalve reserve/ structure relative molecular biomass

% release of gamtes - spawning events
kR = 0.95;
Rspawn = 1;

% *************************************************************************
% ODE 
% *************************************************************************

% options
refine = 6;
options = odeset('Events', @MatRepr,...
    'Maxstep', 0.2, 'OutputSel',1, 'Refine',refine);

%**************************************************************************
% animal
%**************************************************************************

% structure for allowing the integration of several
% other time-dependent variables (idea:
% http://www.mathworks.com/matlabcentral/answers/52740)
x0 = [MV0;ME0;MH0;0]; % initial conditions in a vector
tout = tstart; % the loop should start again when finished for one cycle
xout = x0.'; % the results should be transposed for easier readybility
% and for the possibility to use the data

teout = []; % vector for collecting the time points, where events take place
yeout = []; % vector for the values at the time points of the event
ieout = []; % index i of the event function that vanishes

nTime = length(tspan);
x_list = zeros(nTime, length(x0));
x_list(1,:) = x0;

% calculate the energetic parameter before puberty
for iTime = 2:nTime
    [t,x_,te,ye,ie] = ode45(@LifeCycle, tspan(iTime-1:iTime), x0,...
        options, Tact(iTime), Stationphyto(iTime));
    %     disp('iTime')
    %     disp(iTime)
    nt = length(t);
    tout = [tout; t];
    xout = [xout; x_];
    teout = [teout; te]; % Events at tstart are never reported.
    yeout = [yeout; ye];
    ieout = [ieout; ie];
    if (t(nt) >= tend)
        break;
    end
    x0 = x_(end,:);
    x_list(iTime, :) = x0;
    nt = length(t);
    
    if ie == 1
        juvtend = t(end,:);
        save juv juvtend; % save the day when maturity is reached
    end
    
    
    if ie == 2
        x0 = [x_(nt,1);x_(nt,2);x_(nt,3);0];
        ContReprEventDay = t(end,:);
        % add here a line for collecting the values
        save ContReprEvent ContReprEventDay;      
        
    else
        x0 = x_(end,:);
        x_list(iTime, :) = x0;
    end
end

MV_glob = x_list(:,1); % Bivalve structure biomass
ME_glob = x_list(:,2); % Bivalve reserve biomass
MH_glob = x_list(:,3); % Bivalve maturity investment
MR_glob = x_list(:,4); % Bivalve reproduction buffer

timeJuv = 1:500;
MH_globJuv = MH_glob(1:500);


figure(3)
if MR_glob == 0
    plot(timeJuv, MH_globJuv,'c')
else
    plot(timeJuv, MH_globJuv,'c',teout(1,1), yeout(1,3),'ro')
    h = legend('maturity investment', 'puberty');
    ylabel('biomass in mol C')
    xlabel('time in days')
    xlim([0 500])

    tendnew = tend-1;
    % time for figure 2
    tnewg = 1:tendnew;
    MV_globAd = MV_glob(1:(length(MV_glob)-1));
    ME_globAd = ME_glob(1:(length(ME_glob)-1));
    MR_globAd = MR_glob(1:(length(MR_glob)-1));
    
    
    figure(4)
    plot(tnewg, MV_globAd,'r', tnewg, ME_globAd,'k',tnewg,...
        MR_globAd,'g',teout,yeout(:,4),'ro')
    h = legend('structure biomass', 'reserve biomass','reproduction buffer');
    axis( [0,Inf, 0, Inf] )
    ylabel('biomass in mol C')
    xlabel('time in days')

    % organism lenght
    V = (MV_glob./(dv/wv)).^(1/3); % volumetric length
    uM = 0.297; % shape coefficient
    L = V./uM;
    
    Lg = L(1:(length(L)-1));
    
    figure(5)
    plot(tnewg,Lg);
    ylabel('length in cm')
    xlabel('time in days')
    
  %***Saving*****************************************************************

MV_globMyt = MV_globAd;
ME_globMyt = ME_globAd;
MH_globMyt = MH_globJuv;
MR_globMyt = MR_globAd;
LMyt = Lg;

contReprmatfile = 'contRepr.mat';
save(contReprmatfile,'tnewg', 'MV_globMyt', 'ME_globMyt', 'MR_globMyt', 'MH_globMyt')


% save
contEventRmatfile = 'contEventR.mat';
contReventsMyt = yeout;
contrMR_glob_events = contReventsMyt(:,4); % extract only those numbers of MR where the event occurs
contrJspawnEROut = kR.*contrMR_glob_events./Rspawn; % spawning in in mol C^E C d^-1
contrNspawnOut = contrJspawnEROut/ME0; % calculate the number of gamets
contrTeout = teout;
save(contEventRmatfile, 'contrTeout', 'contReventsMyt', 'contrNspawnOut')

contLengthmatfile = 'contLength.mat';
save(contLengthmatfile,'tnewg', 'LMyt')

load('contRepr.mat','MV_globMyt', 'ME_globMyt', 'MR_globMyt', 'MH_globMyt')
load('contLength.mat','LMyt')
load('contEventR.mat', 'contReventsMyt')

% save the data as csv file
 ResultsResp = [globCRControlTime, globCRControl];
 csvwrite('ResultTableRespContr.dat',ResultsResp);
 type ResultTableRespContr.dat;
 
%**************************************************************************
whos
end
end

function [dX, CR, JxiF, JxiI, Nspawn] = LifeCycle(t, x, Tact, Stationphyto)
% program for evaluating the rates of change of different pysiological
% processes in the study organism

global globCRControlTime
global globCRControl
global globJxiI

globCRControlTime =[globCRControlTime;t];

% ****************************************************
% DEB papameters needed
% ****************************************************
reftempk1 = 20; % reference temperature at which k1 was measured, for all processes
Tref = reftempk1+273; % reference temperature for temperature correction (when TA>0, translate to Kelvin)
T1 = 293; % reference temperature
TA = 5800; % data van der Veer et al. 2006 %5800; % Arrhenius temperature
TL = 275; % 273; % lower temperature boundary range
TH = 296; % 305; % upper temperature boundary range
TAL = 45430; % Arrhenius temperatures for that rate of decrease at both boundaries, data van der Veer et al. 2006
TAH = 31376; % Arrhenius temperatures for that rate of decrease at both boundaries, data van der Veer et al. 2006
ST = (1+exp((TAL/Tact)-(TAL/TL))+ exp((TAH/TH)-(TAH/Tact)))^-1;
ST1 = (1+exp((TAL/T1)-(TAL/TL))+exp((TAH/TH)-(TAH/T1)))^-1;


if Tact >= 305
    Tfact = exp((TA/T1)-(TA/Tact))*ST/ST1; % equation by Kooijmann 2010
elseif Tact < 0
    Tfact = exp((TA/T1)-(TA/Tact))*ST/ST1;
else
    Tfact = exp(TA/Tref - TA/Tact);
end


% values, which were previously saved in ODE file
nxiN = 0.150943396; % mol N / mol C calculated from Redfield ratio % 0.15; % data Saraiva et al 11a % algal nitrogen:carbon ratio in mol N^-1 C
nxiP = 0.009433962; % mol P / mol C % calculated from Redfield ratio %  %0.017 Sakshaug et al 1983 %0.0094; %0.025; % data from graph in Saraiva et al 2012 %0.0094 ;% redfield ratio, check again! % algal phosphorous:carbon ratio in P^-1 C
nEP = 0.006; %Saraiva et al 2012, check again % phosphorous:carbon ratio of bivalve reserve/ structure
nEN = 0.18; % Saraiva et al 2012 %0.238; % Both et al 2012; 0.18; % Saraiva 2012 %0.2381; % Both et al. 2012; 0.22; % data Saraiva et al. 2011a (modelling feeding) %0.18; % nitrogen:carbon ratio of bivalve reserve/ structure
%
%**************************************************************************
% parameters from literature data for Mytilus edulis
%**************************************************************************
YEXV = 0.75; % Yield coefficient of reserves in algal structure in mol C^E

dv = 0.2; % g dw/ cm^3; Rosland et al 09 in Saraiva et al 2011 dv = dE,
          % bivalve structure and reserve specific density
fE = 0.5;% 0.8;  % reserve fraction in algal mass, Saraiva et al 2012
wv = 25.22; % f dw/mol, bivalve reserve/ structure relative molecular biomass

CRm = 0.096; %0.144; % corrected based on data by Wang and Fischer 1996 %0.096; %based on Saravia 2011 clearance rate in m^-3 d^-1 cm^-2 Saraiva et al. (2011b)
JxiFm = 4.8e-4; % algal maximum surface area specific filtration rate
pxiI = 0.9; % algal binding probability, Wang and Fischer, Saraiva et al 2011a
JxiIm = 1.3e4;%1.3e4*500000;%1e112;%0.00355289656902389; %add my pet data  ;%1.3e4; % algal maximum ingestion rate in mol C d^-1, Saraiva et al 2011a
v = 0.056; % energy conductance in cm d^-1, Saraiva et al 2011a
kap = 0.67; % allocation fraction to growth and somatic maintenance, Saraiva et al 2011a
pM = 11.6; % volume specific somatic maintenance in J d^-1 cm-3, Saraiva et al 2011a
EG = 5993; % specific costs for structure in J cm^3, Saraiva et al 2011a

Rspawn = 1; % spawning period in days Saraiva et al 2012

Tspawn = 282.6; % 9.6ï¿½C + 273; % minimum temperatire for spawning in ï¿½C, Hummel et al 1989

uE = 6.97e5; % bivalve reserve chemical potential in J mol^-1, van der Veer et al 2006

kR = 0.95; % reproduction efficiency (Kooijman 2010, Saraiva et al 2012)

kM = (pM/EG); % somatic maintenance rate coefficient d^-1 here = 0.001935591523444
kJ = kM; % volume specific maturity maintenance rate coefficient
MPH = 0.000226686; % cumulative maturity at puberty calculted with Saraiva paper
ME0 = 1.48e-10; % initial reserve mass at optimal food conditions in mol C^E
MVst = dv/wv; % Mstr is the volume-specific structural mass in mol C^V cm ^-3,

%**************************************************************************
% states, which need to be calculated
%**************************************************************************
% The initial conditions are related to the state vector X, where for each
% cumulative stress level over time and in dependency of the differential
% equations values are calculated (state varaiables).
% The vector X contains variables such as maturity, growth, reserve biomass
% etc. (see below)
% In the ODE solver data for certain events are saved, the ODE solver is
% stopped at that occasion and started again with the current values as
% initial conditions. This way the program runs through the different life
% stages of the organism

MV  = x(1); % scaled bivalve structure biomass

%**************************************************************************
% further animal specific parameters that can be calculated
%**************************************************************************

V = (MV/ MVst); % the volume of the organism

YVE = (MVst*uE)/ EG; % Yield coefficient of structure on reserves

%**************************************************************************
% consideration of species specific Arrhenius temperature
% Temperature correction
%**************************************************************************
% k1 is the metabolic rate measured at a certain temperature (reference temperature),
% data should be collected for each metabolic process, replace then k1
% with the corresponding process and add a list of the data here
% included in the model the same reference temperature should be used,
% prepare data with the same equation if just data at other reference
% temperatures are available, Arrhenius temperature is the same for the
% different processes (see details page 18 Kooijman 2010 book

% time of birth assumed to be in May (because water temp is likely to be
% warm enough. However spawning can until Sept and in April water temp
% could also be warm enough already)
% water temperature data from homepage http://www.wassertemperatur.org/sylt


%********************************************************************JEAE*********************************************
%
% Filtration

CR = (CRm/(1+(Stationphyto*CRm/JxiFm)+25.56*CRm/3.5))* (V^(2/3))*Tfact;
globCRControl = [globCRControl; CR];
 
JxiF = CR*Stationphyto*Tfact; % filtration rate in mol C d^-1 g d^-1, CR is scalar (dependent on size)
%


% Ingestion

JxiI = (pxiI*JxiF)/(1+((pxiI*JxiF)/JxiIm))*Tfact;% ingestion rate in mol C d^-1 g d^-1
globJxiI = [globJxiI; JxiI];

%
% Assimilation
JEAV = YEXV*JxiI*(1-fE)*Tfact;%*250;% algal structure assimilation rate in mol C^E d^-1
rC = JxiI*fE; % mol C^E d^-1
rN = JxiI*fE*(nxiN/nEN); % mol C^E d^-1
rP = JxiI*fE*(nxiP/nEP); % mol C^E d^-1

% Algal reserve assimilation rate

JEAE = ((1/rC + 1/rN + 1/rP - (1/(rC+rN))- (1/(rC+rP)) - (1/(rN+rP))+ (1/(rC+rN+rP)))^-1)*Tfact;
JEA = (JEAE + JEAV); % assimilation rate in mol C^E d^-1

% represents the assimilation rate at the reference temperature
%
% faeces production rate not needed for now but maybe later on when
% considering species interactions
% dJpiI = dJxiI - dJEA; % faeces production rate in mol C d ^-1
%
% Somatic maintenance
JES = (pM/uE)*V; % somatice maintenance in mol C^E C d^-1

ME = x(2); % get data for reserve biomass for the corresponding time-step
% Mobilisation
Em = ME/V; % maximum energy density, here the energy density is dependent on size
JEC = (Em/((EG/uE) + kap*Em)*((EG/uE)*v*V^(2/3)+JES)); % Mobilisation flux in mol C^E C d^-1

% bivalve reserve biomass
JE = JEA - JEC;

MH = x(3); % get data for maturity investment for the corresponding time-step

% Growth
JEG = kap*JEC-JES; % flux allocated to growth in mol C^E C d^-1

JVG = YVE*JEG; % Growth in C^V d^-1

if JVG < 0
    JVG = 0;
end

JEJ = kJ*MH; % maturity maintenance, JMER is a state (see above)

% Maturity and reproduction
JER = ((1-kap)*JEC-JEJ); % The flux allocated to maturity and reproduction

if MH < MPH % if MH < MPH, % if maturity is smaller than maturity at puberty, all flux
    % to maturity and reproduction is invested in maturity
    JMER = JER;% flux allocated to reproduction/ maturity in mol C^E C d^-1
    JRER = 0;
else
    JMER = 0;% if maturity at puberty is reached, no flux is invested into
    % the development of maturity
    JRER = JER;% flux allocated to reproduction/ maturity in mol C^E C d^-1
end

MR = x(4);

if Tact >= Tspawn
    JspawnER = kR*MR/ Rspawn; % spawning in in mol C^E C d^-1
    Nspawn = JspawnER/ME0; % number of gamets released in n d^-1
else JspawnER = 0;
    Nspawn = 0;
end

save spawning Nspawn


dX= [JVG; JE; JMER; JRER]; 
end

function [value,isterminal,direction] = MatRepr(t,x_,Tact, Stationphyto)

MH = x_(3);
% This subfunction catches events such as puberty, spawning etc.
% -1 detect zero crossings in the negative direction only
% 0 detect all zero crossings
% 1 detect zero crossings in the positive direction only
% The length of value, isterminal, and direction is the same as the number
% of event functions. The ith element of each vector, corresponds to the
% ith event function.

% data for if cases
% maturity
% MPH = Ehp/uE (Ehp = 1.58e2 J; uE = 6.97e5 J mol^-1) Saraiva et al 2012
MPH = 0.000226686; % cumulative maturity at puberty calculted with Saraiva paper

% spawning

GSRspawn = 0.2; % gonado-somatic ratio to spawn in mol C^R mol^-1 C, Saraiva et al 2012
Tspawn = 282.6; % 9.6ï¿½C + 173; % minimum temperatire for spawning in ï¿½C, Hummel et al 1989

%         %  puberty
%         value(1) = MPH - MH ; % if maturity reaches the level of maturity common at
%         % puberty, the organism starts to invest the energy for maturity and
%         % reproduction (JER) exclusively into reproduction.
%         isterminal(1) = 0; %  stop (1) or not stop the solver at the event
%         % here the solver should be stopped and get provided with new initial
%         % conditions
%         direction(1) = 1; % catch all zero crossings: First the actual maturity is
%         % smaller than maturity at puberty (negative value), then the value
%         % gets positive during puberty. the command direction = 1 instructs
%         % MATLAB to only accept events for which x(2) (that is, x?) is positive
value(1) = x_(3) - MPH; % if the maturity level exceed the normal maturity
% level at puberty, it is assumed that puberty level is reached
isterminal(1) = 1; % stop the solver at the event
direction(1) = 0; % catch all zero crossings when value 1 gets greater
% or smaller than 0#

if Tspawn < Tact
    TspawnAct = 0;
else 
    TspawnAct = 1;
end

% if Tact >= Tspawn
% if the gonado specific ratio (GSR) gets greater than the GSR common at
% spawning, the event spawning is expected to happen
value(2) = (x_(4)/(x_(1)+x_(2)+x_(4)))- GSRspawn + TspawnAct; %Nspawnout - 1; % 2 is just used to ensure a change of
% sign to detect the event
% spwaning event
isterminal(2) = 1; %  stop (1) or not stop the solver at the event
direction(2) = 0; % catch all zero crossings when value 1 gets greater than 0
%  else value(2) = 0; % if itï¿½s too cold spawning does not take place
% unregarding any other circumstances
end



