% (system of differential equations for basic physiological processes and
% processes over the life span of an individual)
% outputs

function MultipleStressorsAddEff

clear
close all
clc

global globCRadd
global globCRaddTime

global globCRControl
global globCRControlTime
global globTact

% initialize these variables to catch the data in a vector
globCRadd = [];
globCRaddTime = [];
globTact = [];

disp('DEB model with additive effects')

% plot(MV_glob,ME_glob)
%%% here add later switch on and off functions
SwitchEnv = 1; % switch for environmental data simulation
% to switch on use 1, to switch off use 0
format compact
tend = 1932; % time in days
tstart = 1; % start at hatching of one juvenile at spring
tspan = [tstart:tend]; % this timespan is used for all calculations

% timedependent changes in the environment alter physiological rates
% first example: temperature
% later on more varables should be added here such as algal concentration
% etc.


%==========================================================================
% temperature
%==========================================================================

StationTemp = 'Temp.xlsx';
sheettemp = 1; % sheet containing the needed temperature values
x1Range_date = 'D2:D1933';
Stationtemp_date = xlsread(StationTemp, sheettemp, x1Range_date);
% delete cells without any data
Stationtemp_date(any(isnan(Stationtemp_date),2),:) = []; % first remove all
% values which have not been measured and contain an empty field, which
% can can create error messages later on
Stationtemp_date = Stationtemp_date + datenum('30DEC1899');

x1Range_values = 'E2:E1933'; %start at 28 April when Temp is > 9,6
Stationtemp_values = xlsread(StationTemp, sheettemp, x1Range_values);
Stationtemp_values(any(isnan(Stationtemp_values),2),:) = []; % first remove all
% values which have not been measured and contain an empty field, which
% can can create error messages later on

% craete file for temperature
Tempmatfile = 'Temp.mat'; % create a file for the stressor
save(Tempmatfile,'Stationtemp_date', 'Stationtemp_values')

load('Temp.mat','Stationtemp_values')
Tact = mat2dataset(Stationtemp_values);
%load('Temp_in_water.mat','Temp_values')
%Tact = mat2dataset(Temp_values);
Tact = double(Tact);
Tact = Tact+273;


% startDate = datenum('21-02-2005');
% endDate = datenum('21-11-2005');
% xData = linspace(startDate, endDate, tspan);
% figure
% plot(xData,Tact)
% plot(tspan,Tact)


% phytoplankton concentration
%%% uptated at 30.08.2017

Station_phyto = 'Phyto2010_NLWKN.xlsx';
sheetphy = 1;  % sheet containing phytoplankton information
x1Range_date = 'C2:C297';
Stationphyto_date = xlsread(Station_phyto, sheetphy, x1Range_date);

% delete cells without any data
Stationphyto_date(any(isnan(Stationphyto_date),2),:) = []; % first remove all
% values which have not been measured and contain an empty field, which
% can can create error messages later on

% values
x1Range_values = 'A2:A297';
Stationphyto_values = xlsread(Station_phyto, sheetphy, x1Range_values);

% Excel and matlab use another system for dates, therfor in ecxel the dates
% have been transformed to a number (use format cells)
% now we need to correct the data in matlab as matlab counts the days since
% the 1 Jänner 0000 and excel counts the days since the 1 Jänner 1900
% therfore it is neccesary to add the 31. Dez 1899 to the imported data
% with te dunction datenum

Stationphyto_date = Stationphyto_date + datenum('30DEC1899');

%--------------------------------------------------------------------------
% TEST
% test if the date is  transferred correctly from excel to matlab
% anser should be 21-Feb-2005
% datestr(StationCd_date(1))
%--------------------------------------------------------------------------

% shows all dates as a string like 21-Feb-2005
Stationphyto_dateStr = datestr(Stationphyto_date(1:end));

phyto_unit = 'mol C/L';
phytomatfile = 'phyto.mat'; % create a file for the stressor

% create a timeseries
Stationphyto_values(any(isnan(Stationphyto_values),2),:) = []; % first remove all
% values which have not been measured and contain an empty field, which
% can can create error messages later on

%--------------------------------------------------------------------------
% TEST
% length(StationCd_values) % use this value later for uncertainty analysis
% StationCdH2O = timeseries(StationCd_values,tspan);
% -------------------------------------------------------------------------

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
% the first 10 lines
tphytoend = tphytoStart + tendnew;
tnew = tphytoStart:tphytoend;

%%%Interpolation%%%
% interpolation of data (phytoplankton is not measured every day and these data
% need to be generated)
%
%%% pchip: Piecewise Cubic Hermite Interpolating Polynomial (PCHIP),
% pchip interpolates using a piecewise cubic polynomial P(x) with these properties:

% On each subinterval xk?x?xk+1, the polynomial P(x) is a cubic Hermite
% interpolating polynomial for the given data points with specified derivatives
% (slopes) at the interpolation points.

% P(x) interpolates y, that is, P(xj)=yj, and the first derivative dPdx is
% continuous. The second derivative d2Pdx2 is probably not continuous so
% jumps at the xj are possible.

% The cubic interpolant P(x) is shape preserving. The slopes at the xj are
% chosen in such a way that P(x) preserves the shape of the data and
% respects monotonicity. Therefore, on intervals where the data is monotonic,
% so is P(x), and at points where the data has a local extremum, so does P(x).
% see https://de.mathworks.com/help/matlab/ref/pchip.html
%%%%%%%%%%%%%%%%%%%

Stationphytopre = interp1(Stationphyto_date,Stationphyto_values,tnew,'pchip');
Stationphyto = Stationphytopre*1000;% calculate with 1000 to convert from liter to m^3


% figure(1)
% subplot(2,1,1)
% plot(Stationphyto_date,Stationphyto_values,'o')
% title('data')
% xlabel('time in days')
% ylabel('phytoplankton in mol C/L ')
% test
% length(StationCd_values)
% length(tnew)
% tspan = 1:2096;
% length(tspan)
% subplot(2,1,2)

% plot(tspan,Stationphyto)
% title('interpolated data')
% xlabel('time in days')
% ylabel('phytoplankton in mol C/L')
%load('phyto.mat','

% little test, which can be activated
% Tact = spline(timeTemp,Tactst,t);
% disp('Tact')
% disp(Tact)
% time = t;
save(phytomatfile,'Stationphyto_date', 'Stationphyto_values', 'phyto_unit', 'Stationphyto');


%

%==========================================================================
% stressors and their interaction factors
%==========================================================================
disp('Switch: environmental conditions (stressors)')
if SwitchEnv == 0
    disp('model is run without stressors')
else
    disp('model is run with stressors')
end

if SwitchEnv == 1
    [EnvData] = EnvironmentalData; % run the function EnvironmentalData to get all
    % environmental data interpolated to daily data
    StationCd = EnvData(1,:);
    StationCd = StationCd';
    
    StationCu = EnvData(2,:);
    StationCu = StationCu';
    
    StationZn = EnvData(6,:);
    StationZn = StationZn';
    
    StationPb = EnvData(4,:);
    StationPb = StationPb';
    
    StationpH = EnvData(5,:);
    StationpH = StationpH';
    
    TactCelsius = Tact-273;
    
    % definition of stressors
    % short description of the stressors, maybe some of these basic data can
    % be used later for an output or grouping, the ginformation about effects is
    % needed for the additive analysis: stressors affecting the same
    % physiological target are assumed to have a different additive effect
    % than two stressors affecting different targets. Based on thi information
    % it is decided if a multiplicatove model (different targets)or a additive
    % model will be used (same targets)
    
    stressor.hazard(1).symbol = 'Cd';
    stressor.hazard(1).name = 'Cadmium';
    stressor.hazard(1).units = 'ug/L'; % add unit types
    stressor.hazard(1).group = 'heavy metals';
    stressor.hazard(1).effects = {'development' 'filtration' 'energy'...
        'growth' 'survival'};
   
  
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
        [t_,YCd_] = ode45(@StressorCd_ode,tspan(iTime-1:iTime),YCd0, [],...
            StationCd(iTime));
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

    %TEST without this stressor influence
    %     eCd = zeros(size(Tact));
    %     eCdadd = zeros(size(Tact));
    
    %stressor.hazard(1).effectadd = CdInterBilanz;
    
    figure(7)
    plot(tnew,eCd)
    title('Cd')
    xlabel('time in days')
    ylabel('eCd')
    
    %         % activate if needed
    %     figure(9)
    %     plot(tnew,eCdadd)
    %     title('Cd add - interaction with temperature and Cu')
    %     xlabel('time in days')
    %     ylabel('eCdadd')
    
    stressor.hazard(1).effectlevel = eCd;
    
  
    
    %**************************************************************************
    % copper
    
    stressor.hazard(2).symbol = 'copper';
    stressor.hazard(2).name = 'Cu';
    stressor.hazard(2).units = 'ug/L'; % add unit types
    stressor.hazard(2).group = 'heavy metals'; %
    stressor.hazard(2).effects = {'filtration' 'growth' 'reproduction' ...
        'energy' 'survival'};
     
    
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
        [t_,YCu_] = ode45(@StressorCu_ode,tspan(iTime-1:iTime),YCu0, [],...
            StationCu(iTime));
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
    
   
    stressor.hazard(2).effectlevel = eCu;
    
    figure(8)
    plot(tnew,eCu)
    title('Cu')
    xlabel('time in days')
    ylabel('eCu')
    
    %**********************************************************************1
    % Zinc
    
    stressor.hazard(3).symbol = 'Zn';
    stressor.hazard(3).name = 'zinc';
    stressor.hazard(3).units = 'unit'; % add unit types
    stressor.hazard(3).group = 'heavy metals';
    stressor.hazard(3).effects = {'filtration' 'energy' 'growth' ...
        'reproduction' 'survival'};

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
        [t_,YZn_] = ode45(@StressorZn_ode,tspan(iTime-1:iTime),YZn0, [],...
            StationZn(iTime));
        
        YZnout = [YZnout;YZn_];
        tout = [tout;t_];
        YZn0 = YZn_(end,:);
        YZn_list(iTime,:) = YZn0;
    end
    
    intZn = YZn_list(:,1);
    
    %normintZn = intZn/C0Zn; % normalised intensity of the stressor (1=
    %tolerance)
    
    eZn = YZn_list(:,3);
    eZn(eZn<0) = 0;
    
      
    % comment to switch off Zn effects
    % eZn = zeros(size(Tact));
    % eZnadd = zeros(size(Tact));
    
    % interaction effects
    stressor.hazard(3).effectlevel = eZn;
    
    
    figure(9)
    plot(tnew,eZn)
    title('Zn')
    xlabel('time in days')
    ylabel('eZn')
    
    %     % activate if needed
    %     figure(10)
    %     plot(tnew,eZnadd)
    %     title('Zn with interactions with Cd and Cu')
    %     xlabel('time in days')
    %     ylabel('eZnadd')
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
        [t_,YPb_] = ode45(@StressorPb_ode,tspan(iTime-1:iTime),YPb0, [],...
            StationPb(iTime),TactCelsius(iTime));
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
       
    % activate for test without Pb effect
    %ePb = zeros(size(Tact));
    
    % exposure time
    % figure(1)
    % plot(tspan,ePb)
    % xlabel('t')
    % ylabel('effect Pb')
    
    % interaction effects

    % final stressor e
    stressor.hazard(4).effectlevel = ePb;
    
    % figures
    figure(10)
    plot(tnew,ePb)
    title('Pb with Temp interaction')
    xlabel('time in days')
    ylabel('ePb')
    
    %     figure(11)
    %     plot(tnew,ePbadd)
    %     title('Pb')
    %     xlabel('time in days')
    %     ylabel('ePbadd')
    
    
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
            StationpH(iTime));
        
        
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
    
    %pH = zeros(size(Tact));%TEST WITHOUT pH inlfuence
    
    % exposure time
    % figure(1)
    % plot(tspan,ePb)
    % xlabel('t')
    % ylabel('effect Pb')
    
    stressor.hazard(5).effectlevel = epH;
    
    figure(11)
    plot(tnew,epH)
    title('pH')
    xlabel('time in days')
    ylabel('epH')
  
    
    %    plot(StationpH,sigmpH)
    %    plot(tspan,epH)
    
    %**************************************************************************
    % Temperature as a stressor
    
    stressor.hazard(6).symbol = 'Temp';
    stressor.hazard(6).name = 'Temp';
    stressor.hazard(6).units = 'temperature in °C'; % add unit types
    stressor.hazard(6).group = 'physical conditions';
    stressor.hazard(6).effects = {'filtration' 'energy' 'development'...
        'growth' 'reproduction' 'survival'};
    
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
    
    %eTemp = zeros(size(Tact));% TEST WITHOUT THE INFLUENE OF TEMP
    
    stressor.hazard(6).effectlevel = eTemp;
    
    figure(12)
    plot(tnew,eTemp)
    title('Temp')
    xlabel('time in days')
    ylabel('eTemp')
    
    %     % activate if needed
    %     plot(tnew,TactCelsius)
    %     title('Temp')
    %     xlabel('time in days')
    %     ylabel('temperature in °C')
    
    % present the data by simple addition of the effect levels (here without
    % considering the targets
    
    % sum the effect levels
    eaddAll = eCd + eCu + eZn + ePb + epH + eTemp;
    
   eaddAllmatfile = 'effectAdd.mat';
    save(eaddAllmatfile,'eaddAll')
    
    % plot
    figure(13)
    plot(tspan,eaddAll)
    xlabel('time in days')
    ylabel('sum of the netto effect levels')
    
    
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
    
    % calculation of the additive effects
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
    
    % calculation of the additive effects
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
    
    % calculation of the additive effects
    sumEffGrow = sum(GrowEff,2);
    sumEffGrow(sumEffGrow<0) = 0;
    
    % % the effects are summed up and weighted by the interaction factors
    % sumEffGrowadd = sumEffGrow.*sumIntGrow;
    Growth = 'growth.mat'; % create a file for the effect
    save(Growth,'sumEffGrow')
    disp('sumEffGrow')
    %     disp(sumEffGrow)
    %
    
    % collecting all stressors affecting the reproductive system
    target = 'reproduction';
    index = arrayfun(@(s)any(strcmp(s.effects, target)), stressor.hazard);
    ReprEff = ([stressor.hazard(index).effectlevel]);
    % ReptInt = ([stressor.hazard(index).interaction]);
    
    % calculation of the additive effects
    sumEffRepr = sum(ReprEff,2);
    sumEffRepr(sumEffRepr<0) = 0;
    % sumIntRepr = sum(ReptInt,2);
    % % the effects are summed up and weighted by the interaction factors
    % sumEffRepradd = sumEffRepr.*sumIntRepr;
    Reproduction = 'reproduction.mat'; % create a file for the effect
    save(Reproduction,'sumEffRepr')
    disp('sumEffRepr')
    disp(sumEffRepr)
    
    %==========================================================================
    % mortality %% use later for Leslie matrix
    %==========================================================================
    %     % mortality is treated as one effect size like the other physiological
    %     % targets. However, this output is used later on in the Leslie model
    %     target = 'mortality';
    %     index = arrayfun(@(s)any(strcmp(s.effects, target)), stressor.hazard);
    %     MortEff = ([stressor.hazard(index).effectlevel]);
    %     % MortInt = ([stressor.hazard(index).interaction]);
    %
    %     % calculation of the additive effects
    %     sumEffMort = sum(MortEff,2);
    %     % sumIntMort = sum(MortInt,2);
    %     % % the effects are summed up and weighted by the interaction factors
    %     % sumEffMortadd = sumEffRepr.*sumIntRepr;
    %
    %     % save additive mortality in a file
    %     Mortality = 'mortality.mat'; % create a file for the effect
    %     save(Mortality,'sumEffMort')
    %     disp('sumEffMort')
    %     disp(sumEffMort)
    %
    %     RespEffadd = ([stressor.hazard(index).effectadd]);
    %     sumEffFiltadd = sum(RespEffadd,2);
    
end
%==========================================================================
% The organism level - DEB model (Saraiva el al 2012)
%==========================================================================

%==========================================================================
% initial conditions
%==========================================================================

MV0 = 3.3e-09; % structural mass at birth in mol C data van der Veer et al 2006;%8.86388638664122e-11; % JVG0;

ME0 = 1.48e-10; % initial reserve mass at optimal food conditions in mol C^E%2.97e-09; % initial biomass, data van der Veer et al 2006

% JER = flux allocated to reproduction and maturity
JMER0 = 0; % 0.0004539; % data from add my pet % or 0.4290; % initial maturity investment = EHb/uE
% (EHb = maturity at birth, uE = bivalve reserve chemical potential)

MH0 = 4.2898e-11; % cumulative maturity at birth, calculated with data 
                  % from Saraiva et al. 2012  Ehb/uE
JRER0 = 0; % 0.0006; % bivalve reproduction buffer
MR0 = JRER0;

MER0 = 0;

%
nxiN = 0.1509; % calculated from Redfield ratio % 0.15; % data Saraiva et al 11a % algal nitrogen:carbon ratio in mol N^-1 C
nxiP = 0.025; % 0.0094 calculated from Redfield ratio % %0.025; % data from graph in Saraiva et al 2012 %0.0094 ;% redfield ratio, check again! % algal phosphorous:carbon ratio in P^-1 C
nEP = 0.006; %Saraiva et al 11b, check again % phosphorous:carbon ratio of bivalve reserve/ structure
nEN = 0.2381; % Both et al. 20120.18; 0.22; % data Saraiva et al. 2011a (modelling feeding) %0.18; % nitrogen:carbon ratio of bivalve reserve/ structure
%
%**************************************************************************
% parameters from literature data
%**************************************************************************
YEXV = 0.75; % Yield coefficient of reserves in algal structure in mol C^E
EHp = 1.58e2; % maturity at puberty in J, Saraiva et al 2011a
dv = 0.2; % g dw/ cm^3; Rosland et al 09 in Saraiva et al 2011 dv = dE, bivalve structure and reserve specific density
fE = 0.5;  % reserve fraction in algal mass, Saraiva et al 2012
pAm = 147.6; % value van der Veer 2006 % 37.99; % maximum specific assimilation rate
TA = 5800; % Arrhenius temperature
reftempk1 = 10; % reference temperature at which k1 was measured, for all processes
Tref = reftempk1+273; % reference temperature for temperature correction (when TA>0, translate to Kelvin)
T1 = 293; % reference temperature
TA = 5800; % data van der Veer et al. 2006 %5800; % Arrhenius temperature
TL = 273; % lower temperature boundary range
TH = 296; % upper temperature boundary range
TAL = 45430; % Arrhenius temperatures for that rate of decrease at both boundaries
TAH = 31376; % Arrhenius temperatures for that rate of decrease at both boundaries
wv = 25.22; % f dw/mol, bivalve reserve/ structure relative molecular biomass
uM = 0.297; % shape coefficient
Rspawn = 1; % Saraiva et al 2012

%
% options
refine = 6;
options = odeset('Events', @MatRepr,...
    'Maxstep', 0.2, 'OutputSel',1, 'Refine',refine);

%
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
%tspan = [tstart:tend]; % timespan
nTime = length(tspan);
x_list = zeros(nTime, length(x0));
x_list(1,:) = x0;
SpawnEventTime = [];

% calculate the energetic parameter before puberty
for iTime = 2:nTime
    [t,x_,te,ye,ie] = ode45(@LifeCycle,tspan(iTime-1:iTime),x0,...
        options,Tact(iTime),Stationphyto(iTime),sumEffFilt(iTime),...
        sumEffGrow(iTime), sumEffRepr(iTime), sumEffEner(iTime));
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


MV_globaddPre = x_list(:,1); % Bivalve structure biomass
ME_globaddPre = x_list(:,2); % Bivalve reserve biomass
MH_globaddPre = x_list(:,3); % Bivalve maturity investment
MR_globaddPre = x_list(:,4); % Bivalve reproduction buffer
% CR_globadd = xout(:,5); % filtration considering additive effects


tnewg = 1:tendnew;
MV_globadd = MV_globaddPre(1:(length(MV_globaddPre)-1));
ME_globadd = ME_globaddPre(1:(length(ME_globaddPre)-1));
MR_globadd = MR_globaddPre(1:(length(MR_globaddPre)-1));

timeJuvadd = 1:500;
MH_globadd = MH_globaddPre(1:500);

%***Saving*****************************************************************
if ie>0
    addReprmatfile = 'AddRepr.mat';
    save(addReprmatfile,'tnewg', 'MV_globAdd', 'ME_globAdd', 'MR_globAdd')
    addEventRmatfile = 'AddEventC.mat';
    addRevents = yeout;
    teoutaddEvents = teout;
    save(addEventRmatfile, 'teoutAddEvents', 'AddRevents')
    
else
    addMatumatfile = 'AddMatu.mat';
    save(addMatumatfile,'tnewg', 'MV_globadd', 'ME_globadd', 'MH_globadd')
end

addReprmatfile = 'addRepr.mat';
save(addReprmatfile,'tnewg', 'MV_globadd', 'ME_globadd', 'MR_globadd', ...
    'MH_globadd')

% release of gamtes - spawning events
kR = 0.95;
Rspawn = 1;
ME0 = 1.48e-10; % initial reserve mass at optimal food conditions in mol C^E
% save
addEventRmatfile = 'AddEventR.mat';
addRevents = yeout;
addMR_glob_events = addRevents(:,4); % extract only those numbers of MR where the event occurs
addJspawnEROut = kR.*addMR_glob_events./Rspawn; % spawning in in mol C^E C d^-1
addNspawnOut = addJspawnEROut/ME0; % calculate the number of gamets
addTeout = teout;

save(addEventRmatfile, 'addTeout', 'addRevents','addNspawnOut')
    V = (MV_globadd./(dv/wv)).^(1/3); % volumetric length
    uM = 0.297; % shape coefficient
    Ladd = V./uM;
    
    % save length data in a matfile
    addLengthmatfile = 'addLength.mat';
    save(addLengthmatfile, 'tnewg', 'Ladd')
    

% do the same for clearance rate data
if length(globCRControl) <= length(globCRadd)
    lengthGraphResp = length(globCRControl);
    globCRaddg = globCRadd(1:lengthGraphResp);
    globCRTimeg = globCRaddTime(1:lengthGraphResp);
    globCRControlg = globCRControl;
else
    lengthGraphResp = length(globCRadd);
    globCRControlg = globCRControl(1:lengthGraphResp);
    globCRTimeg = globCRControlTime(1:lengthGraphResp);
    globCRaddg = globCRadd;
end

% load data from control-scenario
load('contRepr.mat','MV_globMyt', 'ME_globMyt', 'MR_globMyt', 'MH_globMyt')
load('contLength.mat','LMyt')

% ??
AddToMH = NaN*ones(length(MV_globMyt)-length(MH_globMyt),1);
MH_globMytSave = [AddToMH; MH_globMyt];
MH_globaddSave = [AddToMH;MH_globadd];


whos
end

function [dX, CR, JxiF, JxiI, Nspawn, CRpre] = LifeCycle(t, x, Tact, Stationphyto,...
    sumEffFilt, sumEffGrow, sumEffRepr, sumEffEner)
% program for evaluating the rates of change of different pysiological
% processes in the study organism

global globCRadd
global globCRaddTime

globCRaddTime = [globCRaddTime;t];


%TEST
% sumEffFilt = 0;
% sumEffGrow = 0;
% sumEffRepr = 0;
% sumEffEner = 0;

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
nxiN = 0.150943396; % calculated from Redfield ratio % 0.15; % data Saraiva et al 11a % algal nitrogen:carbon ratio in mol N^-1 C
nxiP = 0.009433962; % mol P / mol C % calculated from Redfield ratio % % 0.017; % Sakshaug et al 1983 %0.0094; % calculated from Redfield ratio % %0.025; % data from graph in Saraiva et al 2012 %0.0094 ;% redfield ratio, check again! % algal phosphorous:carbon ratio in P^-1 C
nEP = 0.006; %Saraiva et al 2012, check again % phosphorous:carbon ratio of bivalve reserve/ structure
nEN = 0.18; % Saraiva et al 2012 %0.238; % Both et al 2012; 0.18; % Saraiva 2012 %0.2381; % Both et al. 2012; 0.22; % data Saraiva et al. 2011a (modelling feeding) %0.18; % nitrogen:carbon ratio of bivalve reserve/ structure
%
%**************************************************************************
% parameters from literature data
%**************************************************************************
YEXV = 0.75; % Yield coefficient of reserves in algal structure in mol C^E

EHp = 1.58e2; % maturity at puberty in J, Saraiva et al 2011a
dv = 0.2; % g dw/ cm^3; Rosland et al 09 in Saraiva et al 2011 dv = dE, bivalve structure and reserve specific density
fE = 0.5;  % reserve fraction in algal mass, Saraiva et al 2012
pAm = 147.6; % data by van der Veer 2006 % maximum specific assimilation rate
wv = 25.22; % f dw/mol, bivalve reserve/ structure relative molecular biomass

%**************************************************************************
% parameter Mytilus edulis from Saraiva et al. (2012)
% 'Validation of a Dynamic Energy Budget model
% (DEB) for the blue mussel Mytilus edulis' by Saraiva et al. 2012, Marine
% Ecology Progress Series
%**************************************************************************
CRm = 0.096; % clearance rate in m^-3 d^-1 cm^-2 Saraiva et al. (2011b)
JxiFm = 4.8e-4; % algal maximum surface area specific filtration rate

pxiI = 0.9; % algal binding probability, Saraiva et al 2011a
JxiIm = 1.3e4; % algal maximum ingestion rate in mol C d^-1, Saraiva et al 2011a

v = 0.056; % energy conductance in cm d^-1, Saraiva et al 2011a
kap = 0.67; % allocation fraction to growth and somatic maintenance, Saraiva et al 2011a
pM = 11.6; % volume specific somatic maintenance in J d^-1 cm-3, Saraiva et al 2011a
EG = 5993; % specific costs for structure in J cm^3, Saraiva et al 2011a

Rspawn = 1; % spawning period in days Saraiva et al 2012
GSRspawn = 0.2; % gonado-somatic ratio to spawn in mol C^R mol^-1 C, Saraiva et al 2012
Tspawn = 282.6; % 9.6ï¿½C + 273; % minimum temperatire for spawning in ï¿½C, Hummel et al 1989

uE = 6.97e5; % bivalve reserve chemical potential in J mol^-1, van der Veer et al 2006
kR = 0.95; % reproduction efficiency (Kooijman 2010, Saraiva et al 2012)
Lm = 20; % maximum length

%**************************************************************************
% model based on publication 'Validation of a Dynamic Energy Budget model
% (DEB) for the blue mussel Mytilus edulis' by Saraiva et al. 2012, Marine
% Ecology Progress Series
%**************************************************************************

%**************************************************************************
% states, which need to be calculated
%**************************************************************************
% The initial conditions are related to the state vector X, where for each
% additive stress level over time and in dependency of the differential
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
MPH = 0.000226686; % additive maturity at puberty calculted with Saraiva paper
ME0 = 1.48e-10; % initial reserve mass at optimal food conditions in mol C^E
MVst = dv/wv; % Mstr is the volume-specific structural mass in mol C^V cm ^-3,

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

% % note that a variable that relates to t needs to be directly applied in a
% % differential equation, so it doesnï¿½t work to calculate if t > 30
% % Tenv = 10 etc; and then include this first in Tact = Tenv + 173. The
% % programm can then not relate to t. Therefore this form is needed:

%********************************************************************JEAE*********************************************
%
% Filtration


CRpre = (CRm/(1+((Stationphyto)*CRm/JxiFm)+25.56*CRm/3.5))* (V^(2/3)*Tfact); % Clearance rate in m^3 d^-1, V is scalar
CR = CRpre - CRpre*sumEffFilt;

globCRadd = [globCRadd;CR];

JxiF = CR*Stationphyto*Tfact; % filtration rate in mol C d^-1 g d^-1, CR is scalar (dependent on size)
%

% Ingestion
JxiI = (pxiI*JxiF)/(1+((pxiI*JxiF)/JxiIm))*Tfact; % ingestion rate in mol C d^-1 g d^-1

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

function [value,isterminal,direction] = MatRepr(t,x_,Tact, Stationphyto,...
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
MPH = 0.000226686; % additive maturity at puberty calculted with Saraiva paper

% spawning

GSRspawn = 0.2; % gonado-somatic ratio to spawn in mol C^R mol^-1 C, Saraiva et al 2012
Tspawn = 282.6; % 9.6ï¿½C + 173; % minimum temperatire for spawning in ï¿½C, Hummel et al 1989
kR = 0.95; % reproduction efficiency (Kooijman 2010, Saraiva et al 2012)

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
% unregarding any other ciraddstances
end


