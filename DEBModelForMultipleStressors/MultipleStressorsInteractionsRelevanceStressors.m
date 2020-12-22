
% a program to model the cumulative effects of multiple stressors on blue
% mussels
% (system of differential equations for basic physiological processes and
% processes over the life span of an individual)

clear
close all
clc

format compact
tend = 1932; % time in days
tstart = 1; % start at hatching of one juvenile at spring
tspan = tstart:tend; % this timespan is used for all calculations
stressor = GetStressorIntensities(tspan); 

% run the model with excluding the different stressors one by one to
% determine the relevance of each of the stressors for the modeloutput

% run model without Cd
StressorScene = stressor;
StressorScene.hazard(1).effectlevel = zeros(length(stressor.hazard(1).effectlevel),1); 
[sumEffFilt, sumEffGrow, sumEffRepr, sumEffEner]= RunStressorRelevance(StressorScene);

[dX, CR, JxiF, JxiI, Nspawn, CRpre] = LifeCycle(t, x, Tact, Nneyphyto,...
    sumEffFilt, sumEffGrow, sumEffRepr, sumEffEner);

% compare the different model types (cumulative model and additive model)
% function
% CompareModels

function [stressor] = GetStressorIntensities(tspan)

tend = tspan(end);
% % declare here what you want to investigate
% prompt = 'Select a case compareModels/CompareRelevanceOfStressors:';
% Investigate = input(prompt,'s');

% to compare modeltypes
% Investigate = compareModels;

% to compare the relevance of the stressors
% Investigate = CompareRelevanceOfStressors;

global globCRcum
global globCRcumTime

global globCRControl
global globCRControlTime
global globTact

% initialize these variables to catch the data in a vector
globCRcum = [];
globCRcumTime = [];
globTact = [];

% timedependent changes in the environment alter physiological rates
% first example: temperature

%==========================================================================
% temperature
%==========================================================================

NneyTemp = 'HAMSOM_temp.xlsx';
sheettemp = 1; % sheet containing the needed temperature values
x1Range_date = 'D2:D1933';
Nneytemp_date = xlsread(NneyTemp, sheettemp, x1Range_date);
% delete cells without any data
Nneytemp_date(any(isnan(Nneytemp_date),2),:) = []; % first remove all
% values which have not been measured and contain an empty field, which
% can can create error messages later on
Nneytemp_date = Nneytemp_date + datenum('30DEC1899');

x1Range_values = 'E2:E1933'; %start at 28 April when Temp is > 9,6
Nneytemp_values = xlsread(NneyTemp, sheettemp, x1Range_values);
Nneytemp_values(any(isnan(Nneytemp_values),2),:) = []; % first remove all
% values which have not been measured and contain an empty field, which
% can can create error messages later on

% craete file for temperature
Tempmatfile = 'Temp.mat'; % create a file for the stressor
save(Tempmatfile,'Nneytemp_date', 'Nneytemp_values')

load('Temp.mat','Nneytemp_values')
Tact = mat2dataset(Nneytemp_values);
%load('Temp_in_water.mat','HAMSOMTemp_values')
%Tact = mat2dataset(HAMSOMTemp_values);
Tact = double(Tact);
Tact = Tact+273;

% to correct the dataset (the water temperature cannot be less than zero
% °C)
Tact(Tact<273)=273;

% phytoplankton concentration
Nney_phyto = 'Phyto2010_NLWKN.xlsx';
sheetphy = 1;  % sheet containing phytoplankton information
x1Range_date = 'C2:C297';

%load excel date as number (note the other counting in ecxel than in
%matlab)
Nneyphyto_datePre = xlsread(Nney_phyto, sheetphy, x1Range_date);

% delete cells without any data
Nneyphyto_datePre(any(isnan(Nneyphyto_datePre),2),:) = []; % first remove all
% values which have not been measured and contain an empty field, which
% can can create error messages later on

% Excel and matlab use another system for dates, therfor in ecxel the dates
% have been transformed to a number (use format cells)
% now we need to correct the data in matlab as matlab counts the days since
% the 1st of January 0000 and excel counts the days since the 1st of Jan 1900
% therfore it is neccesary to add the 31. Dez 1899 to the imported data
% with te dunction datenum

Nneyphyto_date = Nneyphyto_datePre + datenum('30DEC1899');

% values
x1Range_values = 'A2:A297';
Nneyphyto_values = xlsread(Nney_phyto, sheetphy, x1Range_values);


% shows all dates as a string like 21-Feb-2005
Nneyphyto_dateStr = datestr(Nneyphyto_date(1:end));

phyto_unit = 'mol C/L';
phytomatfile = 'phyto.mat'; % create a file for the stressor

% create a timeseries
Nneyphyto_values(any(isnan(Nneyphyto_values),2),:) = []; % first remove all
% values which have not been measured and contain an empty field, which
% can can create error messages later on

%**************************************************************************
% interpolate data
%**************************************************************************
% interpolate the data so that they can be used in the other programs
% unit in all programs: days
% matlab method: resample(ts, time, interp_method) (interp_method linear by
% default
tphytoStart = 732430;% corrected date for spawning time as a start for hatching
% (starting time of the model) 28.April 2010 (at 2005 was the temp 9.9
% °C) for Phytoplankton we have only data for one year
tendnew = tend-1; % one needs to be subtracted to get the same timespan as defined above in
% rthe firts 10 lines
tphytoend = tphytoStart + tendnew;
tnew = tphytoStart:tphytoend;

%%%Interpolation%%%
% interpolation of data (phytoplankton is not measured every day and these data
% need to be generated)
%
%%% pchip: Piecewise Cubic Hermite Interpolating Polynomial (PCHIP),
% pchip interpolates using a piecewise cubic polynomial P(x) with these properties:

% On each subinterval, the polynomial P(x) is a cubic Hermite
% interpolating polynomial for the given data points with specified derivatives
% (slopes) at the interpolation points.

% P(x) interpolates y, that is, P(xj)=yj, and the first derivative dPdx is
% continuous. The second derivative d2Pdx2 is probably not continuous so
% jumps at the xj are possible.

% The cubic interpolant P(x) is shape preserving. The slopes at the xj are
% chosen in such a way that P(x) preserves the shape of the data and
% respects monotonicity. Therefore, on intervals where the data is monotonic,
% so is P(x), and at points where the data has a local extremum, so does P(x).
% info copied from https://de.mathworks.com/help/matlab/ref/pchip.html
%%%%%%%%%%%%%%%%%%%

Nneyphytopre = interp1(Nneyphyto_date,Nneyphyto_values,tnew,'pchip');
Nneyphyto = Nneyphytopre*1000;% calculate with 1000 to convert from liter to m^3

%save(phytomatfile,'Nneyphyto_date', 'Nneyphyto_values', 'phyto_unit', 'Nneyphyto')

%%% TEST
% plot(tspan, Nneyphyto)
%tnew = datestr(tnew);
% %
% startDate = datenum('06-01-2005');
% endDate = datenum('11-18-2010');
% xData = linspace(startDate,endDate,1997);
% figure
%plot(xData,Nneyphyto)
%set(gca,'Xtick',xData)
% datetick('x','yy', 'keepsticks');
%endDate =
%plot(tnew,Nneyphyto)

% %TEST
% %shows all dates as a string like 21-Feb-2005
% tnewStr = datestr(tnew);
% disp(tnewStr)

% figure(1)
% subplot(2,1,1)
% plot(Nneyphyto_date,Nneyphyto_values,'o')
% title('data')
% xlabel('time in days')
% ylabel('phytoplankton in mol C/L ')
% test
% length(NneyW1Cd_values)
% length(tnew)
% tspan = 1:2096;
% length(tspan)
% subplot(2,1,2)

% plot(tspan,Nneyphyto)
% title('interpolated data')
% xlabel('time in days')
% ylabel('phytoplankton in mol C/L')
%load('phyto.mat','

% little test, which can be activated
% Tact = spline(timeTemp,Tactst,t);
% disp('Tact')
% disp(Tact)
% time = t;
save(phytomatfile,'Nneyphyto_date', 'Nneyphyto_values', 'phyto_unit', 'Nneyphyto');


% get data from the control-scenario, so that they can be compared to the
% cumulative scenario
disp('DEB model for control-scenario')
SaraivaOdeCTempM(tspan, Tact, Nneyphyto)

load('contRepr.mat','MV_globMyt', 'ME_globMyt', 'MR_globMyt', 'MH_globMyt')
load('contLength.mat','LMyt')


disp('DEB model with cumulative effects')

%==========================================================================
% stressors and their interaction factors
%==========================================================================
% test
% disp('Switch: environmental conditions (stressors)')
% if SwitchEnv == 0
%     disp('model is run without stressors')
% else
%     disp('model is run with stressors')
% end
% %========================================================================
% if SwitchEnv == 1
[EnvData] = EnvironmentalData; % run the function EnvironmentalData to get all
% environmental data interpolated to daily data
NneyW1Cd = EnvData(1,:);
NneyW1Cd = NneyW1Cd';

NneyW1Cu = EnvData(2,:);
NneyW1Cu = NneyW1Cu';

NneyW1O2 = EnvData(3,:);
NneyW1O2 = NneyW1O2';

NneyW1Zn = EnvData(6,:);
NneyW1Zn = NneyW1Zn';

NneyW1Pb = EnvData(4,:);
NneyW1Pb = NneyW1Pb';

NneyW1pH = EnvData(5,:);
NneyW1pH = NneyW1pH';

TactCelsius = Tact-273;

% definition of stressors
% short description of the stressors, maybe some of these basic data can
% be used later for an output or grouping, the ginformation about effects is
% needed for the cumulative analysis: stressors affecting the same
% physiological target are assumed to have a different cumulative effect
% than two stressors affecting different targets. Based on thi information
% it is decided if a multiplicatove model (different targets)or a additive
% model will be used (same targets)

stressor.hazard(1).symbol = 'Cd';
stressor.hazard(1).name = 'Cadmium';
stressor.hazard(1).units = 'ug/L'; % add unit types
stressor.hazard(1).group = 'heavy metals';
stressor.hazard(1).effects = {'development' 'filtration' 'energy'...
    'growth' 'survival'};


%%%%%%%%%%%%%%%%%%%
% Interactions Cd
%%%%%%%%%%%%%%%%%%%
% input for matrix cell:
CdCuInteract = 0.214; % influence of Cu on Cd effect (uptake)(Elliot et al. 1986)
CdZnInteract = 0.252; % influence of Zn on Cd effect (uptake)(Elliot et al. 1986)

% interaction factor derived based on data by George et al. 1983.
CdpHInteract = (-0.1476*NneyW1pH+1.2212);

% based on the publication by Fisher 1986 a formula was generated with
% the cftool to reflect the influence of temperature on the Cd effect
% (here the condition index was observed)
CdTempInteract = (2.644e-05*TactCelsius.^4)-(0.001708*TactCelsius.^3) + (0.04012*TactCelsius.^2)-(0.3779*TactCelsius) + 1.23;

% plot all interactions and save them in a file

CdCuInteractForFigure = ones(length(tnew),1)*CdCuInteract;
CdZnInteractForFigure = ones(length(tnew),1)*CdZnInteract;

figure(8)
plot(tnew,CdpHInteract,'k',tnew,CdTempInteract,'r',tnew,CdCuInteractForFigure,'m',tnew,CdZnInteractForFigure,'b')
datetick('x','mmm yy','keepticks')
xlabel('month and year')
ylabel('CdpH(black),CdTemp(red),CdCu(magenta),CdZn(blue)')
title('Cd_interactions')

newFolder = fullfile(pwd, 'resultsEenv');
if ~exist(newFolder, 'dir')
    mkdir(newFolder)
end

saveas(gcf, fullfile(newFolder, ['CdInteractions', '.fig']));
saveas(gcf, fullfile(newFolder, ['CdInteractionsJPG-', '.jpg']));

disp('CdTempInteract')
disp(CdTempInteract)

% basic data for the initial condition for the differential equations
% initial values
intCd0 = 0; % internal concentration of Cd at reference conditions or at
% the start of the experiment
N0 = 100; % population size or percent funtioning of a process under control conditions
eCd0 = 0; % mortality caused by the stressor
aCd0 = 0; % acclimation to the stressor
YCd0 = [intCd0;N0;eCd0;aCd0];

% here, a program is called, which decribes the effect of a single stressor
% in dependence of exposure time, acclimation, latency time, tolerance
% concentraion etc.
nTime = length(tspan); % use a time, which is suitable for the dataset! Use
% the same time here as in the other ODE solvers
YCd_list = zeros(nTime, length(YCd0));
YCd_list(1, :) = YCd0;

YCdout = []; % vector for collecting the results
tout = []; % vector for collecting the results

for iTime = 2:nTime
    % solving the differential equation
    [t_,YCd_] = ode45(@SchadstoffCd_ode,tspan(iTime-1:iTime),YCd0, [],...
        NneyW1Cd(iTime));
    YCdout = [YCdout;YCd_];
    tout = [tout;t_];
    YCd0 = YCd_(end,:);
    YCd_list(iTime,:) = YCd0;
end

intCd = YCd_list(:,1); % internal concentration after a certain time

% normintCd = intCd/C0Cd; % normalised intensity of the stressor (1= tolerance
% level)
NCd = YCd_list(:,2);
eCd = YCd_list(:,3); % eCd describe the change of effect of Cd over time
eCd(eCd<0) = 0;



%**************************************************************************
% copper

stressor.hazard(2).symbol = 'copper';
stressor.hazard(2).name = 'Cu';
stressor.hazard(2).units = 'ug/L'; % add unit types
stressor.hazard(2).group = 'heavy metals'; %
stressor.hazard(2).effects = {'filtration' 'growth' 'reproduction' ...
    'energy' 'survival'};

%%%%%%%%%%%%%%%%%%%
% Interactions Cu
%%%%%%%%%%%%%%%%%%%
% input for matrix cell:
CuCdInteract = -0.179; % influence of Cd on Cu effect (uptake) (Elliot et al. 1986)
CuZnInteract = 0.698; % influence of Zn  on Cu effect (uptake) (Elliot et al. 1986)

% interaction derived based on data by Akberali et al 1985
exponentpH =-2.674*-NneyW1pH;
CupHInteract = 1.027./(1+6.306e-9*exp(exponentpH));

% interaction derived, based on data by Mubiana et al. 2007
CuTempInteract = 22.45*exp(-0.1606*TactCelsius);

figure(9)
CuZnInteractForFigure = ones(length(tnew),1)*CuZnInteract;
CuCdInteractForFigure = ones(length(tnew),1)*CuCdInteract;

plot(tnew,CupHInteract,'k',tnew,CuTempInteract,'r',tnew,CuCdInteractForFigure,'g',tnew,CuZnInteractForFigure,'b')
datetick('x','mmm yy','keepticks')
xlabel('month and year')
ylabel('interaction factors')
title('Cu_interactions')

newFolder = fullfile(pwd, 'resultsEenv');
if ~exist(newFolder, 'dir')
    mkdir(newFolder)
end

saveas(gcf, fullfile(newFolder, ['CuInteractions', '.fig']));
saveas(gcf, fullfile(newFolder, ['CuInteractionsJPG-', '.jpg']));

figure(11)
CuZnInteractForFigure = ones(length(tnew),1)*CuZnInteract;
CuCdInteractForFigure = ones(length(tnew),1)*CuCdInteract;
plot(tnew,CupHInteract,'k',tnew,CuCdInteractForFigure,'g',tnew,CuZnInteractForFigure,'b')
ylim([-0.3 0.8])

datetick('x','mmm yy','keepticks')
xlabel('month and year')
ylabel('interaction factors')
title('Cu interactions without temp')

newFolder = fullfile(pwd, 'resultsEenv');
if ~exist(newFolder, 'dir')
    mkdir(newFolder)
end

saveas(gcf, fullfile(newFolder, ['CuInteractionsWithoutTemp', '.fig']));
saveas(gcf, fullfile(newFolder, ['CuInteractionsWithoutTempJPG-', '.jpg']));

% initial values
intCu0 = 0;
N0 = 100;
eCu0 = 0;
aCu0 = 0;
YCu0 = [intCu0;N0;eCu0;aCu0];

YCu_list = zeros(nTime, length(YCu0));
YCu_list(1, :) = YCu0;

YCuout = []; % vector for collecting the results
tout = []; % vector for collecting the results

for iTime = 2:nTime
    
    % solving the differential equation
    [t_,YCu_] = ode45(@SchadstoffCu_ode,tspan(iTime-1:iTime),YCu0, [],...
        NneyW1Cu(iTime));
    YCuout = [YCuout;YCu_];
    tout = [tout;t_];
    YCu0 = YCu_(end,:);
    YCu_list(iTime,:) = YCu0;
end

intCu = YCu_list(:,1);

%normintCu = intCu/C0Cu;  % normalised intensity of the stressor (1= tolerance
% level)

eCu = YCu_list(:,3);
eCu(eCu<0) = 0;

% TEST without the influence of Cu
% eCu = zeros(size(Tact));

% figure(2)
% plot(tspan,eCu)
% xlabel('t')
% ylabel('effect Cu')


%**********************************************************************1
% Zinc

stressor.hazard(3).symbol = 'Zn';
stressor.hazard(3).name = 'zinc';
stressor.hazard(3).units = 'unit'; % add unit types
stressor.hazard(3).group = 'heavy metals';
stressor.hazard(3).effects = {'filtration' 'energy' 'growth' ...
    'reproduction' 'survival'};
% stressor.hazard(3).properties = {'hydrophob' 'lipophil' 'persisistent'};


%%%%%%%%%%%%%%%%%%%
% Interactions Zn
%%%%%%%%%%%%%%%%%%%
% input for matrix cell:
ZnCdInteract = 0.28; % influence of Cd on Zn effect (Vercauteren and Blust 1999, mean)
ZnCuInteract = 1.229; % influence of Cu on Zn effect (uptake) (Elliot et al. 1986)

% interaction derived based on data by George 1983
ZnpHInteract = 0.947./(1+3.02e-05*(exp(-1.771*-NneyW1pH)));

% interaction derived based on data by Cotter et al. 1982
ZnTempInteract = 0.001262*TactCelsius.^1.925;


% initial values
intZn0 = 0;
N0 = 100;
eZn0 = 0;
aZn0 = 0;
YZn0 = [intZn0;N0;eZn0;aZn0];

% here, a program is called, which decribes the effect of a single stressor
% in dependence of exposure time, acclimation, latency time, tolerance
% concentraion etc.
nTime = length(tspan); % use a time, which is suitable for the dataset! Use
% the same time here as in the other ODE solvers
YZn_list = zeros(nTime, length(YZn0));
YZn_list(1, :) = YZn0;

YZnout = []; % vector for collecting the results
tout = []; % vector for collecting the results

for iTime = 2:nTime
    % solving the differential equation
    [t_,YZn_] = ode45(@SchadstoffZn_ode,tspan(iTime-1:iTime),YZn0, [],...
        NneyW1Zn(iTime));
    
    YZnout = [YZnout;YZn_];
    tout = [tout;t_];
    YZn0 = YZn_(end,:);
    YZn_list(iTime,:) = YZn0;
end

intZn = YZn_list(:,1);

%normintZn = intZn/C0Zn; % normalised intensity of the stressor 

eZn = YZn_list(:,3);
eZn(eZn<0) = 0;

%
%     %**************************************************************************
% Lead

stressor.hazard(4).symbol = 'Pb';
stressor.hazard(4).name = 'lead';
stressor.hazard(4).units = 'unit'; % add unit types
stressor.hazard(4).group = 'heavy metals';
stressor.hazard(4).effects = {'energy' 'growth' 'reproduction'...
    'survival'}; %
% stressor.hazard(3).properties = {'hydrophob' 'lipophil' 'persisistent'};

%%%%%%%%%%%%%%%%%%%
% Interactions Pb
%%%%%%%%%%%%%%%%%%%
% input for matrix cell:

PbCdInteract = -0.432; % interaction factor derived based on data by Sheir et al 2013

% interaction factor derived based on data by Han Zhao Xiang et al. 2013
PbpHInteract = -0.4506*NneyW1pH+3.752;

% interaction factor derived based on data by Mubiana et al 2007
PbTempInteract = 0.1881*TactCelsius+0.0001;

% initial values
intPb0 = 0;
N0 = 100;
ePb0 = 0;
aPb0 = 0;
YPb0 = [intPb0;N0;ePb0;aPb0];

% here, a program is called, which decribes the effect of a single stressor
% in dependence of exposure time, acclimation, latency time, tolerance
% concentraion etc.
nTime = length(tspan); % use a time, which is suitable for the dataset! Use
% the same time here as in the other ODE solvers
YPb_list = zeros(nTime, length(YPb0));
YPb_list(1, :) = YPb0;

YPbout = []; % vector for collecting the results
tout = []; % vector for collecting the results

% solving the differential equation
for iTime = 2:nTime
    [t_,YPb_] = ode45(@SchadstoffPb_ode,tspan(iTime-1:iTime),YPb0, [],...
        NneyW1Pb(iTime),TactCelsius(iTime));
    YPbout = [YPbout;YPb_];
    tout = [tout;t_];
    YPb0 = YPb_(end,:);
    YPb_list(iTime,:) = YPb0;
end

intPb = YPb_list(:,1);

% normintPb = intPb/C0Pb; % normalised intensity of the stressor (1= tolerance
% level)
%NCd = YCd_list(:,2);
ePb = YPb_list(:,3); % eCd describe the change of effect of Cd with regard to
ePb(ePb<0) = 0;

% integration of the cumulative effects of cumulative effects
% interaction effects


%**************************************************************************
% pH


stressor.hazard(5).symbol = 'pH';
stressor.hazard(5).name = 'pH';
stressor.hazard(5).units = 'no unit'; % add unit types
stressor.hazard(5).group = 'physical conditions';
stressor.hazard(5).effects = {'filtration' 'energy' 'development'...
    'growth' 'reproduction' 'survival'};

% initial values
N0 = 100;
epH0 = 0;
apH0 = 0;
YpH0 = [N0;epH0;apH0];

% here, a program is called, which decribes the effect of a single stressor
% in dependence of exposure time, acclimation, latency time, tolerance
% concentraion etc.
nTime = length(tspan); % use a time, which is suitable for the dataset! Use
% the same time here as in the other ODE solvers
YpH_list = zeros(nTime, length(YpH0));
YpH_list(1, :) = YpH0;

YpHout = []; % vector for collecting the results
tout = []; % vector for collecting the results

% solving the differential equation
for iTime = 2:nTime
    [t_,YpH_] = ode45(@StressorpH_ode,tspan(iTime-1:iTime),YpH0, [],...
        NneyW1pH(iTime));
    
    
    YpHout = [YpHout;YpH_];
    tout = [tout;t_];
    YpH0 = YpH_(end,:);
    YpH_list(iTime,:) = YpH0;
end

% normintPb = intPb/C0Pb; % normalised intensity of the stressor (1= tolerance
% level)
%NCd = YCd_list(:,2);
epH = YpH_list(:,2); % e describes the change of effect of the stressor
epH(epH<0) = 0;

%pH = zeros(size(Tact));% test without pH influence


%**************************************************************************
% Temperature as a stressor

stressor.hazard(6).symbol = 'Temp';
stressor.hazard(6).name = 'Temp';
stressor.hazard(6).units = 'temperature in °C'; % add unit types
stressor.hazard(6).group = 'physical conditions';
stressor.hazard(6).effects = {'filtration' 'energy' 'development'...
    'growth' 'reproduction' 'survival'};

%%%%%%%%%%%%%%%%%%%
% Interactions Temp
%%%%%%%%%%%%%%%%%%%

% interaction derived based on data by Cotter et al 1982
TempZnInteract = 8*10^-5*NneyW1Zn.^2-0.0442*NneyW1Zn+17.053;

% initial values
N0 = 100;
eTemp0 = 0;
aTemp0 = 0;
YTemp0 = [N0;eTemp0;aTemp0];

% here, a program is called, which decribes the effect of a single stressor
% in dependence of exposure time, acclimation, latency time, tolerance
% concentraion etc.
nTime = length(tspan); % use a time, which is suitable for the dataset! Use
% the same time here as in the other ODE solvers
YTemp_list = zeros(nTime, length(YTemp0));
YTemp_list(1, :) = YTemp0;

YTempout = []; % vector for collecting the results
tout = []; % vector for collecting the results

% solving the differential equation
for iTime = 2:nTime
    [t_,YTemp_] = ode45(@StressorTemp_ode,tspan(iTime-1:iTime),YTemp0, [],...
        TactCelsius(iTime));
    YTempout = [tout;t_];
    YTemp0 = YTemp_(end,:);
    YTemp_list(iTime,:) = YTemp0;
end

eTemp = YTemp_list(:,2); % e describes the change of effect of the stressor
eTemp(eTemp<0) = 0;

% integration of the cumulative effects of cumulative effects
% interaction effects
eTempCum = eTemp+(eTemp.*TempZnInteract);

%eTemp = zeros(size(Tact));% TEST WITHOUT THE INFLUENE OF TEMP


%-------------------------------------------------------------------------
% combine interactions
%-------------------------------------------------------------------------
% combine the interactions with consideration of potential double countings
% (e.g. there an increased Cd concentration will alter the effect of Cu
% and an increased Cu concentration on the other hand will alter the effect
% of Cd. However, if the interaction effects would just be added up, the
% effect could be overestimated. Therefore, the effect of both is first
% calculated seperately and then at each time point the sum of both is
% calculated. In a last step, the sum is divided by two. It might be that
% the overlap of the interaction is not evenly distributed. However, as we
% do not know if Cd has a greater or smaller influence on Cu than Cu on Cd,
% we assume an even distribution.

% potential double counting for Cd-Cu interaction

% calculate the sum of both effects for
% Cd ane Cu
CuANDCdInteract = (CdCuInteract + CuCdInteract)./2;
% Cd and Zn
ZnANDCdInteract = (CdZnInteract + ZnCdInteract)./2;
% Cu and Zn
CuANDZnInteract = (CuZnInteract + ZnCuInteract)./2;
% Zn and Temp
ZnANDTempInteract = (ZnTempInteract + TempZnInteract)./2;


% integration of the cumulative effects of cumulative effects

%***************
% Cd
%***************
eCdCum = eCd+(eCd.*CdCuInteract)+(eCd.*CdZnInteract)+(eCd.*CdpHInteract)+(eCd.*CdTempInteract);

%TEST without this stressor influence
%     eCd = zeros(size(Tact));
%     eCdCum = zeros(size(Tact));

%stressor.hazard(1).effectCum = CdInterBilanz;

figure(8)
plot(tnew,eCd,'k',tnew,eCdCum,'r')
datetick('x','mmm yy','keepticks')
xlabel('month and year')
ylabel('eCd and eCdCum')
title('Cd')

newFolder = fullfile(pwd, 'resultsEenv');
if ~exist(newFolder, 'dir')
    mkdir(newFolder)
end

saveas(gcf, fullfile(newFolder, ['CdECum', '.fig']));
saveas(gcf, fullfile(newFolder, ['CdECumJPG-', '.jpg']));

%         % activate if needed
%     figure(9)
%     plot(tnew,eCdCum)
%     title('Cd cum - interaction with temperature and Cu')
%     xlabel('time in days')
%     ylabel('eCdCum')

eCdCum(eCdCum<0) = 0;
stressor.hazard(1).effectlevel = eCdCum;


%***************
% Cu
%***************
eCuCum = eCu+(eCu.*CuCdInteract)+(eCu.*CuZnInteract)+(eCu.*CupHInteract)+(eCu.*CuTempInteract);

stressor.hazard(2).effectlevel = eCuCum;

figure(9)
plot(tnew,eCu,'k',tnew,eCuCum,'r')
datetick('x','mmm yy','keepticks')
title('Cu')
xlabel('month and year')
ylabel('eCu')

saveas(gcf, fullfile(newFolder, ['CuECum', '.fig']));
saveas(gcf, fullfile(newFolder, ['CuECumJPG-', '.jpg']));

%*************************
% Zn
%*************************
% interaction effects

eZnCum = eZn+(eZn.*ZnCdInteract)+(eZn.*CuZnInteract)+(eZn.*ZnpHInteract)+ (eZn.*ZnTempInteract);

% comment to swith off Zn effects
% eZn = zeros(size(Tact));
% eZnCum = zeros(size(Tact));

% interaction effects
stressor.hazard(3).effectlevel = eZnCum;


figure(10)
plot(tnew,eZn,'k',tnew,eZnCum,'r')
datetick('x','mmm yy','keepticks')
title('Zn')
xlabel('month and year')
ylabel('eZn')

saveas(gcf, fullfile(newFolder, ['ZnECum', '.fig']));
saveas(gcf, fullfile(newFolder, ['ZnECumJPG-', '.jpg']));

%     % activate if needed
%     figure(10)
%     plot(tnew,eZnCum)
%     title('Zn with interactions with Cd and Cu')
%     xlabel('time in days')
%     ylabel('eZnCum')
%***************
% Pb
%***************

ePbCum = ePb+(ePb.*PbCdInteract)+(ePb.*PbpHInteract)+(ePb.*PbTempInteract);

ePbCum(ePbCum<0) = 0;

% activate for test without Pb effect
%ePb = zeros(size(Tact));

% exposure time
% figure(1)
% plot(tspan,ePb)
% xlabel('t')
% ylabel('effect Pb')

% interaction effects
% interaction with Temperature is included in the uptake process!!

% final stressor e
stressor.hazard(4).effectlevel = ePbCum;

% figures
figure(11)
plot(tnew,ePb,'k',tnew,ePbCum,'r')
datetick('x','mmm yy','keepticks')
title('Pb')
xlabel('month and year')
ylabel('ePb')

saveas(gcf, fullfile(newFolder, ['PbECum', '.fig']));
saveas(gcf, fullfile(newFolder, ['PbECumJPG-', '.jpg']));

%     figure(11)
%     plot(tnew,ePbCum)
%     title('Pb')
%     xlabel('time in days')
%     ylabel('ePbCum')



%********************
% pH
%********************
stressor.hazard(5).effectlevel = epH;

figure(12)
plot(tnew,epH,'k')
datetick('x','mmm yy','keepticks')
title('pH')
xlabel('month and year')
ylabel('epH')

saveas(gcf, fullfile(newFolder, ['pHE', '.fig']));
saveas(gcf, fullfile(newFolder, ['pHEJPG-', '.jpg']));

%******************
% Temp
%******************

stressor.hazard(6).effectlevel = eTempCum;

figure(13)
plot(tnew,eTemp,'k',tnew,eTempCum,'r')
datetick('x','mmm yy','keepticks')
title('Temp')
xlabel('month and year')
ylabel('eTemp')

saveas(gcf, fullfile(newFolder, ['TempE', '.fig']));
saveas(gcf, fullfile(newFolder, ['TempEJPG-', '.jpg']));


%     % activate if needed
%     plot(tnew,TactCelsius)
%     title('Temp')
%     xlabel('time in days')
%     ylabel('temperature in °C')

% present the data by simple addition of the effect levels (here without
% considering the targets

%==========================================================================
% function CalcCumulativeEffects[]= (stressor.hazard.effectlevel)
%-----------------------------------------------------------------
% sum the effect levels
eCumAll = eCdCum + eCuCum + eZnCum + ePbCum + epH + eTempCum;


% plot
figure(14)
plot(tspan,eCumAll,'r', tspan, eaddAll,'k')
datetick('x','mmm yy','keepticks')
xlabel('month and year')
ylabel('sum of the netto effect levels')

%save the current figure in a new folder
newFolder = fullfile(pwd, 'results');
if ~exist(newFolder, 'dir')
    mkdir(newFolder)
end

saveas(gcf, fullfile(newFolder, ['ComparisonAddCum', '.fig']));
saveas(gcf, fullfile(newFolder, ['ComparisonAddCum-', '.jpg']));


end
%--------------------------------------------------------------------------
 
function CompareModels

    % run the model for additive effects
    disp('DEB model with additive effects')
    SaraivaOdeAddEff
    
    load('addRepr.mat','MV_globadd', 'ME_globadd', 'MR_globadd','MH_globadd')
    load('addLength.mat','Ladd')
    
    % load data
    % all effects summed up without interaction effects
    load('effectAdd.mat','eaddAll')
end


function [sumEffFilt, sumEffGrow, sumEffRepr, sumEffEner]= RunStressorRelevance(stressor)
%**************************************************************************
% TARGETS
%**************************************************************************
% here the stressors are sorted according to the effect they cause
% in mybiosis there will be a structure for sorting the different effects
% the same target does not refer to a molecular target. Instead, the method
% is transferred to the higher physiological level of physiological systems
% This is done for simplifying the concept because the concept should be
% applicable to all kinds of organisms and it should be possible to
% structure the literature information in a transparant way


% collecting all stressors affecting the filtration system of the organism
target = 'filtration';
index = arrayfun(@(s)any(strcmp(s.effects, target)), stressor.hazard);
FiltEff = ([stressor.hazard(index).effectlevel]);
%RespInt = ([stressor.hazard(index).interaction]);

% calculation of the cumulative effects
sumEffFilt = sum(FiltEff,2);
sumEffFilt(sumEffFilt<0) = 0;
filtration = 'filtration.mat'; % create a file for the effect
save(filtration,'sumEffFilt')
disp('sumEffFilt')
disp(sumEffFilt)


% collecting all stressors affecting the endocrine system of the organism
target = 'energy';
index = arrayfun(@(s)any(strcmp(s.effects, target)), stressor.hazard);
EnerEff = ([stressor.hazard(index).effectlevel]);
%EnerInt = ([stressor.hazard(index).interaction]);

% calculation of the cumulative effects
sumEffEner = sum(EnerEff,2);
sumEffEner(sumEffEner<0) = 0;
Energy = 'energy.mat'; % create a file for the effect
save(Energy,'sumEffEner')
disp('sumEffEner')
disp(sumEffEner)


% collecting all stressors affecting the growth
target = 'growth';
index = arrayfun(@(s)any(strcmp(s.effects, target)), stressor.hazard);
GrowEff = ([stressor.hazard(index).effectlevel]);
%GrowInt = ([stressor.hazard(index).interaction]);

% calculation of the cumulative effects
sumEffGrow = sum(GrowEff,2);
sumEffGrow(sumEffGrow<0) = 0;

% sumIntGrow = sum(GrowInt,2);
% % the effects are summed up and weighted by the interaction factors
% sumEffGrowCum = sumEffGrow.*sumIntGrow;
Growth = 'growth.mat'; % create a file for the effect
save(Growth,'sumEffGrow')
disp('sumEffGrow')
%     disp(sumEffGrow)
%
%
% collecting all stressors affecting the reproductive system
target = 'reproduction';
index = arrayfun(@(s)any(strcmp(s.effects, target)), stressor.hazard);
ReprEff = ([stressor.hazard(index).effectlevel]);
% ReptInt = ([stressor.hazard(index).interaction]);

% calculation of the cumulative effects
sumEffRepr = sum(ReprEff,2);
sumEffRepr(sumEffRepr<0) = 0;
% sumIntRepr = sum(ReptInt,2);
% % the effects are summed up and weighted by the interaction factors
% sumEffReprCum = sumEffRepr.*sumIntRepr;
Reproduction = 'reproduction.mat'; % create a file for the effect
save(Reproduction,'sumEffRepr')
disp('sumEffRepr')
disp(sumEffRepr)

%==========================================================================
% mortality %% use later for development of Leslie matrix
%==========================================================================
%     % mortality is treated as one effect size like the other physiological
%     % targets. However, this output is used later on in the Leslie model
%     target = 'mortality';
%     index = arrayfun(@(s)any(strcmp(s.effects, target)), stressor.hazard);
%     MortEff = ([stressor.hazard(index).effectlevel]);
%     % MortInt = ([stressor.hazard(index).interaction]);
%
%     % calculation of the cumulative effects
%     sumEffMort = sum(MortEff,2);
%     % sumIntMort = sum(MortInt,2);
%     % % the effects are summed up and weighted by the interaction factors
%     % sumEffMortCum = sumEffRepr.*sumIntRepr;
%
%     % save cumulative mortality in a file
%     Mortality = 'mortality.mat'; % create a file for the effect
%     save(Mortality,'sumEffMort')
%     disp('sumEffMort')
%     disp(sumEffMort)
%
%     RespEffCum = ([stressor.hazard(index).effectCum]);
%     sumEffFiltCum = sum(RespEffCum,2);


%==========================================================================
% The organism level - DEB model (Saraiva el al 2012)
%==========================================================================

%==========================================================================
% initial conditions
%==========================================================================

MV0 = 3.3e-09; % structural mass at birth in mol C data van der Veer et al 2006;%8.86388638664122e-11; % JVG0;
ME0 = 1.48e-10; % initial reserve mass at optimal food conditions in mol C^E%2.97e-09; % initial biomass, data van der Veer et al 2006
MH0 = 4.2898e-11; % cumulative maturity at birth, calculated with data
% from Saraiva et al. 2012  Ehb/uE

%**************************************************************************
% parameters from literature data
%**************************************************************************

dv = 0.2; % g dw/ cm^3; Rosland et al 09 in Saraiva et al 2011 dv = dE, bivalve structure and reserve specific density

wv = 25.22; % f dw/mol, bivalve reserve/ structure relative molecular biomass

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
xout = x0.'; % the results should be transposed for easier readability
% and for the possibility to use the data

teout = []; % vector for collecting the time points, where events take place
yeout = []; % vector for the values at the time points of the event
ieout = []; % index i of the event function that vanishes

nTime = length(tspan);
x_list = zeros(nTime, length(x0));
x_list(1,:) = x0;
SpawnEventTime = [];

% calculate the energetic parameter before puberty
for iTime = 2:nTime
    [t,x_,te,ye,ie] = ode45(@SaraivaT,tspan(iTime-1:iTime),x0,...
        options,Tact(iTime),Nneyphyto(iTime),sumEffFilt(iTime),...
        sumEffGrow(iTime), sumEffRepr(iTime), sumEffEner(iTime));
    %     disp('iTime')
    %     disp(iTime)
    nt = length(t);
    tout = [tout; t];
    xout = [xout; x_];
    teout = [teout; te]; % events at tstart are never reported.
    yeout = [yeout; ye];
    ieout = [ieout; ie];
    if (t(nt) >= tend)
        break;
    end
    x0 = x_(end,:);
    x_list(iTime, :) = x0;
    nt = length(t);
    
    if ie == 1
        % save the day of the maturation event
        juvtend = t(end,:);
        save juv juvtend;
    end
    
    
    if ie == 2
        x0 = [x_(nt,1);x_(nt,2);x_(nt,3);0];
    else
        x0 = x_(end,:);
        x_list(iTime, :) = x0;
    end
end


MV_globCumPre = x_list(:,1); % Bivalve structure biomass
ME_globCumPre = x_list(:,2); % Bivalve reserve biomass
MH_globCumPre = x_list(:,3); % Bivalve maturity investment
MR_globCumPre = x_list(:,4); % Bivalve reproduction buffer
% CR_globCum = xout(:,5); % filtration considering cumulative effects


tnewg = 1:tendnew;
MV_globCum = MV_globCumPre(1:(length(MV_globCumPre)-1));
ME_globCum = ME_globCumPre(1:(length(ME_globCumPre)-1));
MR_globCum = MR_globCumPre(1:(length(MR_globCumPre)-1));


timeJuvCum = 1:500;
MH_globCum = MH_globCumPre(1:500);

%***Saving*****************************************************************
if ie>0
    cumReprmatfile = 'cumRepr.mat';
    save(cumReprmatfile,'tnewg', 'MV_globCum', 'ME_globCum', 'MR_globCum')
    cumEventRmatfile = 'cumEventC.mat';
    cumRevents = yeout;
    teoutCumEvents = teout;
    save(cumEventRmatfile, 'teoutCumEvents', 'cumRevents')
    
else
    cumMatumatfile = 'cumMatu.mat';
    save(cumMatumatfile,'tnewg', 'MV_globCum', 'ME_globCum', 'MH_globCum')
end

cumReprmatfile = 'cumRepr.mat';
save(cumReprmatfile,'tnewg', 'MV_globCum', 'ME_globCum', 'MR_globCum', ...
    'MH_globCum')

% release of gamtes - spawning events
kR = 0.95;
Rspawn = 1;
ME0 = 1.48e-10; % initial reserve mass at optimal food conditions in mol C^E
% save
cumEventRmatfile = 'cumEventR.mat';
cumRevents = yeout;
cumMR_glob_events = cumRevents(:,4); % extract only those numbers of MR where the event occurs
cumJspawnEROut = kR.*cumMR_glob_events./Rspawn; % spawning in in mol C^E C d^-1
cumNspawnOut = cumJspawnEROut/ME0; % calculate the number of gamets
cumTeout = teout;

save(cumEventRmatfile, 'cumTeout', 'cumRevents','cumNspawnOut')

%**************************************************************************
% figures
%**************************************************************************
disp('figures showing the result for the cumulative model alone')
figure(4)
if MR_globCum == 0
    plot(timeJuvCum, MH_globCum,'c')
else
    plot(timeJuvCum, MH_globCum,'c', teout(1,1), yeout(1,3),'ro')
    h = legend('maturity investment', 'puberty');
    title('szenario: cumulative pressures')
    ylabel('biomass in mol C')
    xlabel('time in days')
    xlim([0 500])

    %save the current figure in a new folder
    newFolder = fullfile(pwd, 'results');
    if ~exist(newFolder, 'dir')
        mkdir(newFolder)
    end
    saveas(gcf, fullfile(newFolder, ['SummaryEaryDevelopment', '.fig']));
    saveas(gcf, fullfile(newFolder, ['SummaryEaryDevelopment-', '.jpg']));
    
    figure(5)
    plot(tnewg, MV_globCum,'r', tnewg, ME_globCum,'k', tnewg, ...
        MR_globCum,'g', teout,yeout(:,4),'ro'),%teout,yeout(:,4),'ro' )
    h = legend('structure biomass', 'reserve biomass', 'reproduction buffer');
    axis( [0,Inf, 0, Inf] )
    title('szenario: cumulative pressures')
    ylabel('biomass in mol C')
    xlabel('time in days')
    
    %save the current figure in a new folder
    newFolder = fullfile(pwd, 'results');
    if ~exist(newFolder, 'dir')
        mkdir(newFolder)
    end
    
    saveas(gcf, fullfile(newFolder, ['SummaryReproduction', '.fig']));
    saveas(gcf, fullfile(newFolder, ['SummaryReproduction-', '.jpg']));
       
    % organism lenght
    V = (MV_globCum./(dv/wv)).^(1/3); % volumetric length
    uM = 0.297; % shape coefficient
    LCum = V./uM;
    
    cumLengthmatfile = 'cumLength.mat';
    save(cumLengthmatfile, 'tnewg', 'LCum')
       
    
    figure(6)
    plot(tnewg,LCum);
    title('growth, szenario: cumulative pressure')
    ylabel('length in cm')
    xlabel('time in days')
    xlim([0 2100])
    ylim([0 8])
    
    %save the current figure in a new folder
    newFolder = fullfile(pwd, 'results');
    if ~exist(newFolder, 'dir')
        mkdir(newFolder)
    end
    
    saveas(gcf, fullfile(newFolder, ['GrowthBlue', '.fig']));
    saveas(gcf, fullfile(newFolder, ['GrowthBlue-', '.jpg']));
    
end


% do the same for clearance rate data
if length(globCRControl) <= length(globCRcum)
    lengthGraphResp = length(globCRControl);
    globCRcumg = globCRcum(1:lengthGraphResp);
    globCRTimeg = globCRcumTime(1:lengthGraphResp);
    globCRControlg = globCRControl;
else
    lengthGraphResp = length(globCRcum);
    globCRControlg = globCRControl(1:lengthGraphResp);
    globCRTimeg = globCRControlTime(1:lengthGraphResp);
    globCRcumg = globCRcum;
end


if Investigate == 'compareModels'
    load('addRepr.mat','MV_globadd', 'ME_globadd', 'MR_globadd','MH_globadd')
    load('addLength.mat','Ladd')
    
    figure(7)
    plot(tnewg, MV_globCum,'r', tnewg, MV_globMyt, 'b', tnewg, MV_globadd, 'k')
    title('Bivalve structure biomass')
    ylabel('biomass in mol C')
    xlabel('time in days')
    
    %save the current figure in a new folder
    newFolder = fullfile(pwd, 'results');
    if ~exist(newFolder, 'dir')
        mkdir(newFolder)
    end
    
    saveas(gcf, fullfile(newFolder, ['StructureBiomass', '.fig']));
    saveas(gcf, fullfile(newFolder, ['StructureBiomass-', '.jpg']));
    
    figure(8)
    plot(tnewg, ME_globCum, 'r',tnewg, ME_globMyt , 'b', tnewg, ME_globadd, 'k')
    title('Bivalve reserve biomass')
    ylabel('biomass in mol C')
    xlabel('time in days')
    
    %save the current figure in a new foldernewFolder = fullfile(pwd, 'results');
    if ~exist(newFolder, 'dir')
        mkdir(newFolder)
    end
    
    saveas(gcf, fullfile(newFolder, ['reserves', '.fig']));
    saveas(gcf, fullfile(newFolder, ['reserves-', '.jpg']));
    
    figure(9)
    tnewJuv = 1:500;
    plot(tnewJuv, MH_globCum, 'r',tnewJuv, MH_globMyt , 'b', tnewJuv, MH_globadd, 'k')
    title('Bivalve maturity investment')
    xlim([0 500])
    ylabel('biomass in mol C')
    xlabel('time in days')
    
    %save the current figure in a new folder
    newFolder = fullfile(pwd, 'results');
    if ~exist(newFolder, 'dir')
        mkdir(newFolder)
    end
    
    saveas(gcf, fullfile(newFolder, ['maturity', '.fig']));
    saveas(gcf, fullfile(newFolder, ['maturity-', '.jpg']));
    
    figure(10)
    plot(tnewg, MR_globCum, 'r',tnewg, MR_globMyt , 'b', tnewg,MR_globadd, 'k')
    title('Bivalve reproduction buffer')
    ylabel('biomass in mol C')
    xlabel('time in days')
    
    %save the current figure in a new folder
    newFolder = fullfile(pwd, 'results');
    if ~exist(newFolder, 'dir')
        mkdir(newFolder)
    end
    
    saveas(gcf, fullfile(newFolder, ['reproduction', '.fig']));
    saveas(gcf, fullfile(newFolder, ['reproduction-', '.jpg']));
    
    figure(11)
    plot(tnewg,LCum,'r',tnewg,LMyt,'b', tnewg,Ladd,'k')
    %plot(tout,Lcum);
    title('growth')
    ylabel('length in cm')
    xlabel('time in days')
    xlim([0 2100])
    ylim([0 10])
    
    %save the current figure in a new folder
    newFolder = fullfile(pwd, 'results');
    if ~exist(newFolder, 'dir')
        mkdir(newFolder)
    end
    
    saveas(gcf, fullfile(newFolder, ['growth', '.fig']));
    saveas(gcf, fullfile(newFolder, ['growth-', '.jpg']));
    
    figure(12)
    plot(globCRTimeg,globCRcumg, 'r',globCRTimeg,globCRControlg, 'b')
    title('clearance rate')
    ylabel('clearance rate in m^3 per day')
    xlabel('time in days')
    
    
    %save the current figure in a new folder
    newFolder = fullfile(pwd, 'results');
    if ~exist(newFolder, 'dir')
        mkdir(newFolder)
    end
    
    saveas(gcf, fullfile(newFolder, ['filtration', '.fig']));
    saveas(gcf, fullfile(newFolder, ['filtration-', '.jpg']));
    
       
    % correct the length of the vector for later calculations
    AddToMH = NaN*ones(length(MV_globMyt)-length(MH_globMyt),1);
    
    MH_globMytSave = [AddToMH; MH_globMyt];
    MH_globCumSave = [AddToMH;MH_globCum];
    MH_globAddSave = [AddToMH;MH_globadd];
    
   
    % ResultMatrixData = [tnewg', MV_globMyt, MV_globCum, ME_globMyt,...
    %     ME_globCum, MH_globMytSave, MH_globCumSave, MR_globMyt, MR_globCum, LMyt, LCum];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate how much % change is due to additive and cumulative effects
    % calculate the sum of all result vectors
    % control-scenario
    SumMVCont = sum(MV_globMyt);
    SumMECont = sum(ME_globMyt);
    SumMHCont = nansum(MH_globMytSave); %NaNs are expected ignored for calculating ths sum
    SumMRCont = nansum(MR_globMyt); %NaNs are expected ignored for calculating ths sum
    SumLCont = sum(LMyt);
    load('contEventR.mat','contrNspawnOut')
    SumNspawnCont = sum(contrNspawnOut);
    
    % additive scenario
    SumMVAdd = sum(MV_globadd);
    SumMEAdd = sum(ME_globadd);
    SumMHAdd = nansum(MH_globAddSave); %NaNs are expected ignored for calculating ths sum
    SumMRAdd = nansum(MR_globadd); %NaNs are expected ignored for calculating ths sum
    SumLAdd = sum(Ladd);
    load('AddEventR.mat','addNspawnOut')
    SumNspawnAdd = sum(addNspawnOut);
    
    % cumulative senario
    SumMVCum = sum(MV_globCum);
    SumMECum = sum(ME_globCum);
    SumMHCum = nansum(MH_globCumSave); %NaNs are expected ignored for calculating ths sum
    SumMRCum = nansum(MR_globCum); %NaNs are expected ignored for calculating ths sum
    SumLCum = sum(LCum);
    load('cumEventR.mat','cumNspawnOut')
    SumNspawnCum = sum(cumNspawnOut);
    
    % effect strength
    SumEadd = sum(eaddAll);
    SumEcum = sum(eCumAll);
    
    % calculate the percent difference for the comparison between control and
    % additive effects
    changeMVAdd = (SumMVAdd*100/SumMVCont)-100;
    changeMEAdd = (SumMEAdd*100/SumMECont)-100;
    changeMHAdd = (SumMHAdd*100/SumMHCont)-100;
    changeMRAdd = (SumMRAdd*100/SumMRCont)-100;
    changeLAdd = (SumLAdd*100/SumLCont)-100;
    changeNspawnAdd = (SumNspawnAdd*100/SumNspawnCont)-100;
    
    % calculate the percent difference for the comparison between control and
    % cumulative effects
    changeMVCum = (SumMVCum*100/SumMVCont)-100;
    changeMECum = (SumMECum*100/SumMECont)-100;
    changeMHCum = (SumMHCum*100/SumMHCont)-100;
    changeMRCum = (SumMRCum*100/SumMRCont)-100;
    changeLCum = (SumLCum*100/SumLCont)-100;
    changeNspawnCum = (SumNspawnCum*100/SumNspawnCont)-100;
    
    % calculate the percent difference for the comparison between additive and
    % cumulative effects
    changeMVCumAdd = (SumMVCum*100/SumMVAdd)-100;
    changeMECumAdd = (SumMECum*100/SumMEAdd)-100;
    changeMHCumAdd = (SumMHCum*100/SumMHAdd)-100;
    changeMRCumAdd = (SumMRCum*100/SumMRAdd)-100;
    changeLCumAdd = (SumLCum*100/SumLAdd)-100;
    changeNspawnCumAdd = (SumNspawnCum*100/SumNspawnAdd)-100;
    
    changeEcumAdd = (SumEcum*100/SumEadd)-100;
    
    
    % save the results data in an excelfile
    ResultMatrixData = [tnewg', MV_globMyt, MV_globadd, MV_globCum,...
        ME_globMyt, ME_globadd, ME_globCum, MH_globMytSave, MH_globAddSave,...
        MH_globCumSave, MR_globMyt, MR_globadd, MR_globCum, LMyt, Ladd,...
        LCum];
    
    % save as csv file
    csvwrite('ResultMatrixDEBMyt.dat',ResultMatrixData)
    type ResultMatrixDEBMyt.dat
    
    % add the headlines
    cHeader = {'time', 'MV_glob','MV_globadd' 'MV_globCum', 'ME_glob',...
        'ME_globadd','ME_globCum', 'MH_glob','MH_globadd', 'MH_globCum',...
        'MR_glob','MR_globadd','MR_globCum', 'LMyt','Ladd','LCum'};
    
    % insert commas
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))];
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); %cHeader in text with commas
    
    % open the csv file and write the header to the file
    % (The 'w' input specifies write access, the 'n' input specifies native byte
    % ordering, and 'Shift_JIS' specifies the character encoding scheme.)
    
    fid= fopen('ResultMatrixDEBMyt.dat','w');
    fprintf(fid,'%s\n',textHeader);
    fclose(fid);
    
    % write data to end of file
    dlmwrite('ResultMatrixDEBMyt.dat',ResultMatrixData,'-append');
    
    % do the same for the data for filtration
    % write a csv file for filtration
    ResultsRespData = [globCRTimeg, globCRControlg, globCRcumg];
    csvwrite('ResultRespDEBMyt.dat',ResultsRespData)
    type ResultRespDEBMyt.dat
    
    % add header
    cHeaderResp = {'time', 'globCRControlg', 'globCRcumg'};
    % insert commas
    commaHeaderResp = [cHeaderResp;repmat({','},1,numel(cHeaderResp))];
    commaHeaderResp = commaHeaderResp(:)';
    textHeaderResp = cell2mat(commaHeaderResp); %cHeader in text with commas
    
    % open the csv file and write the header to the file
    % (The 'w' input specifies write access, the 'n' input specifies native byte
    % ordering, and 'Shift_JIS' specifies the character encoding scheme.)
    
    fid= fopen('ResultRespDEBMyt.dat','w');
    fprintf(fid,'%s\n',textHeaderResp);
    fclose(fid);
    
    % write data to end of file
    dlmwrite('ResultRespDEBMyt.dat',ResultsRespData,'-append');
    
    %******************************************
    
    % results cumulative effects (strength)
    ResultsEcum = [eaddAll,eCumAll];
    %s save as csv file
    csvwrite('ResultsEcum.dat',ResultsEcum)
    type ResultsEcum.dat
    % add the headlines
    cHeader = {'eaddAll', 'eCumAll'};
    
    % insert commas
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))];
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); %cHeader in text with commas
    
    % open the csv file and write the header to the file
    % (The 'w' input specifies write access, the 'n' input specifies native byte
    % ordering, and 'Shift_JIS' specifies the character encoding scheme.)
    
    fid= fopen('ResultsEcum.dat','w');
    fprintf(fid,'%s\n',textHeader);
    fclose(fid);
    
    % write data to end of file
    dlmwrite('ResultsEcum.dat',ResultsEcum,'-append');
    
    %****************************************************************
    % results comparisons between contro, additive and cumulative scenario
    
    ResultsComp = [changeMVAdd,changeMEAdd,changeMHAdd,...
        changeMRAdd,changeLAdd,changeMVCum,changeMECum,changeMHCum,...
        changeMRCum,changeLCum,changeMVCumAdd,changeMECumAdd,...
        changeMHCumAdd,changeMRCumAdd,changeLCumAdd,changeEcumAdd,...
        changeNspawnAdd,changeNspawnCum,changeNspawnCumAdd];
    
    csvwrite('ResultsComp.dat',ResultsComp)
    type ResultsComp.dat
    
    % add the headlines
    cHeader = {'changeMVAdd','changeMEAdd','changeMHAdd',...
        'changeMRAdd','changeLAdd','changeMVCum','changeMECum','changeMHCum',...
        'changeMRCum','changeLCum','changeMVCumAdd','changeMECumAdd',...
        'changeMHCumAdd','changeMRCumAdd','changeLCumAdd','changeEcumAdd',...
        'changeNspawnAdd','changeNspawnCum','changeNspawnCumAdd'};
    
    % insert commas
    commaHeader = [cHeader;repmat({','},1,numel(cHeader))];
    commaHeader = commaHeader(:)';
    textHeader = cell2mat(commaHeader); %cHeader in text with commas
    
    % open the csv file and write the header to the file
    % (The 'w' input specifies write access, the 'n' input specifies native byte
    % ordering, and 'Shift_JIS' specifies the character encoding scheme.)
    
    fid= fopen('ResultsComp.dat','w');
    fprintf(fid,'%s\n',textHeader);
    fclose(fid);
    
    % write data to end of file
    dlmwrite('ResultsComp.dat',ResultsComp,'-append');
end

whos
end

function [dX, CR, JxiF, JxiI, Nspawn, CRpre] = LifeCycle(t, x, Tact, Nneyphyto,...
    sumEffFilt, sumEffGrow, sumEffRepr, sumEffEner)
% program for evaluating the rates of change of different pysiological
% processes in the study organism

global globCRcum
global globCRcumTime

globCRcumTime = [globCRcumTime;t];

%TEST
% sumEffFilt = 0;
% sumEffGrow = 0;
% sumEffRepr = 0;
% sumEffEner = 0;

reftempk1 = 20; % reference temperature at which k1 was measured, for all processes, Saraiva et al. 2012
Tref = reftempk1+273; % reference temperature for temperature correction (when TA>0, translate to Kelvin)
TA = 7022; % data van der Veer et al. 2006 %5800; % Arrhenius temperature
TL = 275; % lower temperature boundary range
TH = 296; %  upper temperature boundary range
TAL = 45430; % Arrhenius temperatures for rate of decrease at lower boundaries, data van der Veer et al. 2006
TAH = 31376; % Arrhenius temperatures for rate of decrease at upper boundaries, data van der Veer et al. 2006
ST = (1+exp((TAL/Tact)-(TAL/TL))+ exp((TAH/TH)-(TAH/Tact)))^-1;
ST1 = (1+exp((TAL/Tref)-(TAL/TL))+exp((TAH/TH)-(TAH/Tref)))^-1;


if Tact >= 305
    Tfact = exp((TA/Tref)-(TA/Tact))*ST/ST1; % equation by Kooijmann 2010
elseif Tact < 0
    Tfact = exp((TA/Tref)-(TA/Tact))*ST/ST1;
else
    Tfact = exp(TA/Tref - TA/Tact);
end


% values, which were previously saved in ODE file
nxiN = 0.150943396; % calculated from Redfield ratio % 0.15; % data Saraiva et al 11a % algal nitrogen:carbon ratio in mol N^-1 C
nxiP = 0.009433962; % mol P / mol C % calculated from Redfield ratio 
nEP = 0.006; %Saraiva et al 2012, phosphorous:carbon ratio of bivalve reserve/ structure
nEN = 0.18; % Saraiva et al 2012 %0.238; % Both et al 2012; 0.18; % Saraiva 2012 %0.2381; % Both et al. 2012; 0.22; % data Saraiva et al. 2011a (modelling feeding) %0.18; % nitrogen:carbon ratio of bivalve reserve/ structure
%
%**************************************************************************
% parameters from literature data
%**************************************************************************
YEXV = 0.75; % Yield coefficient of reserves in algal structure in mol C^E

dv = 0.2; % g dw/ cm^3; Rosland et al 09 in Saraiva et al 2011 dv = dE, bivalve structure and reserve specific density
fE = 0.5;  % reserve fraction in algal mass, Saraiva et al 2012

wv = 25.22; % f dw/mol, bivalve reserve/ structure relative molecular biomass

%**************************************************************************
% initial conditions
%**************************************************************************

%**************************************************************************
% parameter Mytilus edulis from Saraiva et al. (2012)
% 'Validation of a Dynamic Energy Budget model
% (DEB) for the blue mussel Mytilus edulis' by Saraiva et al. 2012, Marine
% Ecology Progress Series
%**************************************************************************
CRm = 0.096; % Maximum surface area specific clearance rate in m^-3 d^-1 cm^-2 Saraiva et al. (2011b)
JxiFm = 4.8e-4; % algal maximum surface area specific filtration rate

pxiI = 0.9; % algal binding probability, changed. in Saraiva et al 2011a value of 0.4

JxiIm = 1.3e4; % algal maximum ingestion rate in mol C d^-1, Saraiva et al 2011a

v = 0.056; % energy conductance in cm d^-1, Saraiva et al 2011a
kap = 0.67; % allocation fraction to growth and somatic maintenance, Saraiva et al 2011a
pM = 11.6; % volume specific somatic maintenance in J d^-1 cm-3, Saraiva et al 2011a
EG = 5993; % specific costs for structure in J cm^3, Saraiva et al 2011a

Rspawn = 1; % spawning period in days Saraiva et al 2012

Tspawn = 282.6; % 9.6ï¿½C + 273; % minimum temperatire for spawning in ï¿½C, Hummel et al 1989

uE = 6.97e5; % bivalve reserve chemical potential in J mol^-1, van der Veer et al 2006

kR = 0.95; % reproduction efficiency (Kooijman 2010, cited in Saraiva et al 2012)


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
% all in all itï¿½s a kind of feedback system

MV  = x(1); % scaled bivalve structure biomass


%**************************************************************************
% animal specific parameters
%**************************************************************************
%
kM = (pM/EG); % somatic maintenance rate coefficient d^-1 here = 0.001935591523444
kJ = kM; % volume specific maturity maintenance rate coefficient
MPH = 0.000226686; % cumulative maturity at puberty calculted with Saraiva paper
ME0 = 1.48e-10; % initial reserve mass at optimal food conditions in mol C^E
MVst = dv/wv; % Mstr is the volume-specific structural mass in mol C^V cm ^-3,

V = (MV/ MVst); % the volume of the organism

YVE = (MVst*uE)/ EG; % Yield coefficient of structure on reserves

% denotes the conversion coefficient from structural volume to C-mole, see
% Kooijman 2010, page 84
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

CRpre = (CRm/(1+((Nneyphyto)*CRm/JxiFm)+25.56*CRm/3.5))* (V^(2/3)*Tfact); % Clearance rate in m^3 d^-1, V is scalar
CR = CRpre - CRpre*sumEffFilt;

globCRcum = [globCRcum;CR];

JxiF = CR*Nneyphyto*Tfact; % filtration rate in mol C d^-1 g d^-1, CR is scalar (dependent on size)
%

% Ingestion
% one = ones(size(L)); % use this trick if L is a vector

JxiI = (pxiI*JxiF)/(1+((pxiI*JxiF)/JxiIm))*Tfact; % ingestion rate in mol C d^-1 g d^-1

%
% Assimilation
JEAV = YEXV*JxiI*(1-fE)*Tfact; % algal structure assimilation rate in mol C^E d^-1
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
Empre = ME/V; % maximum energy density, here the energy density is dependent on size
Em = Empre - Empre*sumEffEner;

JEC = (Em/((EG/uE) + kap*Em)*((EG/uE)*v*V^(2/3)+JES)); % Mobilisation flux in mol C^E C d^-1
%dJEC = min(dJEC,dJEA);
% bivalve reserve biomass
JE = JEA - JEC;

MH = x(3); % get data for maturity investment for the corresponding time-step

% Growth
JEG = kap*JEC-JES; % flux allocated to growth in mol C^E C d^-1

JVGpre = YVE*JEG; % Growth in C^V d^-1
JVG = JVGpre - (sumEffGrow*JVGpre);

if JVG < 0
    JVG = 0;
end

JEJ = kJ*MH; % maturity maintenance, JMER is a state (see above)
% JEJpre = kJ*MH; % maturity maintenance, JMER is a state (see above)
% JEJ = JEJpre - (sumEffDevel*JEJpre);


% Maturity and reproduction
JER = ((1-kap)*JEC-JEJ); % The flux allocated to maturity and reproduction

if MH < MPH % if MH < MPH, % if maturity is smaller than maturity at puberty, all flux
    % to maturity and reproduction is invested in maturity
    JMER = JER;% flux allocated to reproduction/ maturity in mol C^E C d^-1
    JRER = 0;
else
    JMER = 0;% if maturity at puberty is reached, no flux is invested into
    % the development of maturity
    %JRER = JER;
    JRERpre = JER;% flux allocated to reproduction/ maturity in mol C^E C d^-1
    JRER = JRERpre - JRERpre*sumEffRepr;
end


% disp('JMER')
% disp(JMER)
% Spawning

% disp('JRER')
% disp(JRER)
MR = x(4);

if Tact >= Tspawn
    JspawnER = kR*MR/ Rspawn; % spawning in in mol C^E C d^-1
    Nspawn = JspawnER/ME0; % number of gamets released in n d^-1
else JspawnER = 0;
    Nspawn = 0;
end

save spawning Nspawn


dX = [JVG; JE; JMER; JRER];
end

function [value,isterminal,direction] = MatRepr(t,x_,Tact, Nneyphyto,...
    sumEffFilt, sumEffGrow, suumEffRepr, sumEffEner)

global globTact
globTact = [globTact;Tact];

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
MPH = 0.000226686; % cumulative maturity at puberty calculted with Saraiva
% paper (see above), MPH = Ehp/uE

% spawning
% GSR gonado somatic ratio in mol C^R mol^-1 C
% GSR = MR_glob/(MV_glob + ME_glob + MR_glob);
GSRspawn = 0.2; % gonado-somatic ratio to spawn in mol C^R mol^-1 C, Saraiva et al 2012
Tspawn = 282.6; % 9.6ï¿½C + 173; % minimum temperatire for spawning in ï¿½C, Hummel et al 1989
% kR = 0.95; % reproduction efficiency (Kooijman 2010, Saraiva et al 2012)

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


