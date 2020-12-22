function [EnvData] = EnvironmentalData

% delete(findall(0,'Type','figure')) % delete all figures

% Here data from monitoring programs can be entered, sorted, interpolated
% and saved for later use in the DEB model program
% time series data are saved in the format dd-mmm-yyyy HH:MM:SS according
% one of the accepted data-string formats of matlab (seconds are not
% relevant and were not measures, however as the format needs these data we
% simply always used 00 for the seconds. The time of the day might indeed
% be relevant for some of the environmental data)


% define tspan
tendpre = 1931; % 3287; % number of days the model runs
tstart = 732528; % date of possible spawning based on Temperature values 
% of the HAMSOM model (tranlated into a Matlabdate number) 
tend = 732528 + tendpre; % first number represents the date for spawning time as a start for hatching
% (starting time of the model)
tnewStress = tstart:tend;
TimeEnv = 1:length(tnewStress);

%==========================================================================
%O2
%==========================================================================
% Norderney data
%%%% NneyW1_O2 = 'NneyW1.xlsx'; 
NneyW1_O2 = 'NneyW1O2AllYears.xlsx';
sheetO2 = 2; % sheet containing all information concerning O2



% load the relevant data from the excel file
% dates

x1Range_date = 'I2:I188';
NneyW1O2_date = xlsread(NneyW1_O2, sheetO2, x1Range_date);

% delete cells without any data 
NneyW1O2_date(any(isnan(NneyW1O2_date),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

% values
x1Range_values = 'X2:X188';

NneyW1O2_values = xlsread(NneyW1_O2, sheetO2, x1Range_values);
% delete also here cels without data
NneyW1O2_values(any(isnan(NneyW1O2_values),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

% Excel and matlab use another system for dates, therfor in ecxel the dates
% have been transformed to a number (use format cells)
% now we need to correct the data in matlab as matlab counts the days since
% the 1 Jänner 0000 and excel counts the days since the 1 Jänner 1900
% therfore it is neccesary to add the 31. Dez 1899 to the imported data
% with te dunction datenum

NneyW1O2_date = NneyW1O2_date + datenum('30DEC1899');

%--------------------------------------------------------------------------
% TEST
% disp('O2 dates')
% disp(NneyW1O2_date)
% test if the date is  transferred correctly from excel to matlab 
% datestr(NneyW1O2_date(1))
% shows all dates as a string like 21-Feb-2005
% NneyW1O2_dateStr = datestr(NneyW1O2_date(1:end));
%--------------------------------------------------------------------------

O2_in_water_unit = 'mg/L';
O2matfile = 'O2_in_water.mat'; % create a file for the stressor
save(O2matfile,'NneyW1O2_date', 'NneyW1O2_values', 'O2_in_water_unit') 

% create a timeseries 
NneyW1O2_values(any(isnan(NneyW1O2_values),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

%--------------------------------------------------------------------------
% TEST
% length(NneyW1Cd_values) % use this value later for uncertainty analysis
% NneyW1CdH2O = timeseries(NneyW1Cd_values,tspan);
% -------------------------------------------------------------------------

%**************************************************************************
% interpolate data
%**************************************************************************
% interpolate the data so that they can be used in the other programs
% unit in all programs: days
% matlab method: resample(ts, time, interp_method) (interp_method linear by
% default

% TEST
% disp(NneyW1O2_date)
% disp(NneyW1O2_values)

%%NneyW1O2 = interp1(NneyW1O2_date,NneyW1O2_values,tnewStress,'pchip');

% figure(1)
% plot(tnewStress,NneyW1O2,'k',NneyW1O2_date,NneyW1O2_values,'ro')
% %%plot(tnewStress,NneyW1O2)
% title('O2')

% save(O2matfile,'NneyW1O2_date', 'NneyW1O2_values', 'O2_in_water_unit', 'NneyW1O2')
% % save NneyW1O2
% load(O2matfile,'NneyW1O2_date', 'NneyW1O2_values', 'O2_in_water_unit', 'NneyW1O2')

save(O2matfile,'NneyW1O2_date', 'NneyW1O2_values', 'O2_in_water_unit')
% save NneyW1O2
load(O2matfile,'NneyW1O2_date', 'NneyW1O2_values', 'O2_in_water_unit')

% modelling seasonal pattern with the curve fitting tool and a fourier
% model
% NOTE: to make the data available in the cftool one needs to get the data
% through the command window (just copy and paste the corresponding lines
% where the data are derived from excel)
[fitresult, gof] = createFitO2(NneyW1O2_date, NneyW1O2_values);

% modelling with cftool
function [fitresult, gof] = createFitO2(NneyW1O2_date, NneyW1O2_values)
%CREATEFIT(NNEYW1O2_DATE,NNEYW1O2_VALUES)
%  Create a fit.
%
%  Data for 'O2' fit:
%      X Input : NneyW1O2_date
%      Y Output: NneyW1O2_values
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 29-Aug-2019 12:24:42


% Fit: 'O2'.
[xData, yData] = prepareCurveData( NneyW1O2_date, NneyW1O2_values );

% Set up fittype and options.
% simulation of a yearly cycle
% w = 0.0172 because 2pi/365 = 0.0172 (365 days), see method in Fidino and
% Magle 2017 Fourier series
ft = fittype( 'a0+a1*cos(x*0.0172) + b1*sin(x*0.0172)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.964888535199277 0.157613081677548 0.970592781760616];


% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Create a figure for the plots.
figure(1)
title('O2')
%figure( 'Name', 'O2' );

% Plot fit with data.
subplot( 2, 1, 1 );
plot( fitresult, xData, yData, 'predobs' );
% Label axes
dateO2 = datestr(NneyW1O2_date,'mmm/dd');
xlabel dateO2
ylabel 'O2 in mg/L'
ax = gca;
ax.XTick = xData;
datetick('x','mmm yy')
grid on

% Plot residuals.
subplot( 2, 1, 2 );
plot( fitresult, xData, yData, 'residuals' );
% Label axes
xlabel dateO2
ylabel 'O2 in mg/L'
ax = gca;
ax.XTick = xData;
datetick('x','mmm yy')

grid on
end

%  show the results
disp('model for O2')
disp(fitresult)
disp(gof)

resultVec = [];
for x= tnewStress
    resultVecAdd = fitresult(x);
    resultVec = [resultVec,resultVecAdd];
end

NneyW1O2 = resultVec;


% plot the results with regard to the time period chosen (model repeats
% anyway each year, therefore graph not shown explicitly)
% subplot( 3, 1, 3 );
% plot(tnewStress,resultVec)
% dateO2new = datestr(NneyW1O2_date,'mmm/dd');
% xlabel dateO2new
% ylabel 'O2 in mg/L'
% title 'model result'
% ax = gca;
% ax.XTick = xData;
% datetick('x','mmm yy')
% grid on

% -------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% enable figures if needed
% figure(1)
% % subplot(2,1,1)
% plot(NneyW1O2_date,NneyW1O2_values,'o',tnewStress,NneyW1O2,'r')
% title('data and interpolated data O2')
% xlabel('time in days')
% ylabel('O2 in mg/L')
% % test 
% % length(NneyW1Cd_values)
% xlabel('time in days')
% ylabel('CO2 in mg/L')
% -------------------------------------------------------------------------

% *************************************************************************
% create a  time series for later use (e.g. manipulation etc.)
%**************************************************************************
NneyW1O2_values_ts = timeseries(NneyW1O2_values,NneyW1O2_date,'name',...
    'NneyW1O2');

%==========================================================================
%pH
%==========================================================================

tendpre = 1931; % 3287; % number of days the model runs
tstart = 732528; % date of possible spawning based on Temperature values 
% of the HAMSOM model (tranlated into a Matlabdate number) 
tend = 732528 + tendpre; % first number represents the date for spawning time as a start for hatching
% (starting time of the model)
tnewStress = tstart:tend;

% Norderney data
NneyW1_pH = 'NneyW1pHAllYears.xlsx'; 
sheetpH = 3; % sheet containing all information concerning pH

% possible test if Matlab gets the sheet names
% [status,sheets] = xlsfinfo(NneyW1_Cd )

%NneyW1Cd = xlsread(NneyW1_Cd); % load all data from the excel file

% load the relevant data from the excel file
% dates
x1Range_date = 'L2:L937';
NneyW1pH_date = xlsread(NneyW1_pH, sheetpH, x1Range_date);

% delete cells without any data 
NneyW1pH_date(any(isnan(NneyW1pH_date),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

% values
x1Range_values = 'S2:S937';
NneyW1pH_values = xlsread(NneyW1_pH, sheetpH, x1Range_values);

% delete cells without any data 
NneyW1pH_values(any(isnan(NneyW1pH_values),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

% Excel and matlab use another system for dates, therfor in ecxel the dates
% have been transformed to a number (use format cells)
% now we need to correct the data in matlab as matlab counts the days since
% the 1 Jänner 0000 and excel counts the days since the 1 Jänner 1900
% therfore it is neccesary to add the 31. Dez 1899 to the imported data
% with te dunction datenum

NneyW1pH_date = NneyW1pH_date + datenum('30DEC1899');

%--------------------------------------------------------------------------
% TEST
% test if the date is  transferred correctly from excel to matlab 
% anser should be 21-Feb-2005
% datestr(NneyW1Cd_date(1))
%--------------------------------------------------------------------------

% shows all dates as a string like 21-Feb-2005
NneyW1pH_dateStr = datestr(NneyW1pH_date(1:end));


pH_unit = ' '; % Since pH is a logarithmic scale, a difference 
% of one pH unit is equivalent to a tenfold difference in hydrogen ion concentration.

%--------------------------------------------------------------------------
% TEST
% length(NneyW1Cd_values) % use this value later for uncertainty analysis
% NneyW1CdH2O = timeseries(NneyW1Cd_values,tspan);
% -------------------------------------------------------------------------

%**************************************************************************
% interpolate data
%**************************************************************************
% interpolate the data so that they can be used in the other programs
% unit in all programs: days
% matlab method: resample(ts, time, interp_method) (interp_method linear by
% default

% TEST
% disp(NneyW1Cd_date)
% disp(NneyW1Cd_values)

% interpolation
NneyW1pH = interp1(NneyW1pH_date,NneyW1pH_values,tnewStress,'nearest');
pHmatfile = 'pH_.mat'; % create a file for the stressor
save(pHmatfile,'NneyW1pH_date', 'NneyW1pH_values', 'pH_unit', 'NneyW1pH') 

figure(2)
plot(tnewStress,NneyW1pH,'k', NneyW1pH_date, NneyW1pH_values,'ro')
title('pH')
ax = gca;
ax.XTick = tnewStress;
datetick('x','mmm yy')
axis([min(tnewStress) max(tnewStress) 7.5 max(NneyW1pH_values)])


save NneyW1pH
% -------------------------------------------------------------------------
%TEST
%shows all dates as a string like 21-Feb-2005
%tnewStr = datestr(tnew);
%disp(tnewStr)
%figure(1)
%plot(tnew,neyW1Cd_values_interp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enable figures if needed
% figure(5)
% % subplot(2,1,1)
% % subplot(2,1,1)
% plot(NneyW1Cu_date,NneyW1Cu_values,'o',tnewStress,NneyW1Cu,'r')
% title('data and interpolated data pH')
% xlabel('time in days')
% ylabel('pH')
% % test 
% % length(NneyW1Cd_values)

% ------------------------------------------------------------------------

% *************************************************************************
% create a  time series for later use (e.g. manipulation etc.)
%**************************************************************************
NneyW1pH_values_ts = timeseries(NneyW1pH_values,NneyW1pH_date,'name',...
    'NneyW1pH');

% %======================================================
% % data from water column measurements

%==========================================================================
% Cd
%==========================================================================

% Norderney data
NneyW1_Cd = 'NneyW1.xlsx'; 
sheetCd = 1; % sheet containing all information concerning Cd


% possible test if Matlab gets the sheet names
% [status,sheets] = xlsfinfo(NneyW1_Cd )

%NneyW1Cd = xlsread(NneyW1_Cd); % load all data from the excel file

% load the relevant data from the excel file
% dates
x1Range_date = 'H2:H13';
NneyW1Cd_date = xlsread(NneyW1_Cd, sheetCd, x1Range_date);

% delete cells without any data 
NneyW1Cd_date(any(isnan(NneyW1Cd_date),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

% values
x1Range_values = 'P2:P13';

% load the relevant data from excel file
NneyW1Cd_values = xlsread(NneyW1_Cd, sheetCd, x1Range_values);

% delete also here cels without data
NneyW1Cd_values(any(isnan(NneyW1Cd_values),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

%%%%%%%%%%%%%%%%%%TEST%%%%%%%%%%%%%%%%%%%%%


% Excel and matlab use another system for dates, therfor in ecxel the dates
% have been transformed to a number (use format cells)
% now we need to correct the data in matlab as matlab counts the days since
% the 1 Jänner 0000 and excel counts the days since the 1 Jänner 1900
% therfore it is neccesary to add the 31. Dez 1899 to the imported data
% with te dunction datenum

NneyW1Cd_date = NneyW1Cd_date + datenum('30DEC1899');

%--------------------------------------------------------------------------
% TEST
% test if the date is  transferred correctly from excel to matlab 
% datestr(NneyW1Cd_date(1));
%--------------------------------------------------------------------------

% shows all dates as a string like 21-Feb-2005
NneyW1Cd_dateStr = datestr(NneyW1Cd_date(1:end));

Cd_in_water_unit = 'ug/L';

% create a timeseries 
NneyW1Cd_values(any(isnan(NneyW1Cd_values),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

%--------------------------------------------------------------------------
% TEST

% plot(NneyW1Cd_date,NneyW1Cd_values)
% length(NneyW1Cd_values) % use this value later for uncertainty analysis
% NneyW1CdH2O = timeseries(NneyW1Cd_values,tspan);
% -------------------------------------------------------------------------

%**************************************************************************
% interpolate data
%**************************************************************************
% interpolate the data so that they can be used in the other programs
% unit in all programs: days
% matlab method: resample(ts, time, interp_method) (interp_method linear by
% default

% TEST
% disp(NneyW1Cd_date)
% disp(NneyW1Cd_values)

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


NneyW1Cd = interp1(NneyW1Cd_date,NneyW1Cd_values,tnewStress,'pchip');
% disp('NneyW1Cdpre')
% disp(NneyW1Cdpre)
%NneyW1Cdpr = NneyW1Cdpre([100:tend]);
% disp('Cdtime')
% disp(length(tnewStress))
% disp('Cd')
% disp(length(NneyW1Cd))
% plot(tnewStress,NneyW1Cd)
% disp('NneyW1Cdpr')
% disp(NneyW1Cdpr)

Cumatfile = 'Cd_in_water.mat'; % create a file for the stressor
save(Cumatfile,'NneyW1Cd_date', 'NneyW1Cd_values', 'Cd_in_water_unit', 'NneyW1Cd')
% save NneyW1Cd
figure(3)
plot(tnewStress,NneyW1Cd,'k',NneyW1Cd_date,NneyW1Cd_values,'ro')
datetick('x','mmm yy','keepticks')
xlabel('month and year')
ylabel('Cd concentration in \mug/L')
title('Cd')
load('Cd_in_water.mat','NneyW1Cd_date', 'NneyW1Cd_values', 'Cd_in_water_unit', 'NneyW1Cd')

%save the current figure in a new folder
    newFolder = fullfile(pwd, 'resultsEnv');
    if ~exist(newFolder, 'dir')
        mkdir(newFolder)
    end
    
    saveas(gcf, fullfile(newFolder, ['CdConc', '.fig']));
    saveas(gcf, fullfile(newFolder, ['CdConcJPG-', '.jpg']));
% -------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enable figures if needed
% figure(3)
% % subplot(2,1,1)
% plot(NneyW1Cd_date,NneyW1Cd_values,'o',tnewStress,NneyW1Cd,'r')
% title('data and interpolated data Cd')
% xlabel('time in days')
% ylabel('Cd in ug/L')
% % test 
% % length(NneyW1Cd_values)
% xlabel('time in days')
% ylabel('Cd in ug/L')
% -------------------------------------------------------------------------

% *************************************************************************
% create a  time series for later use (e.g. manipulation etc.)
%**************************************************************************
NneyW1Cd_values_ts = timeseries(NneyW1Cd_values,NneyW1Cd_date,'name',...
    'NneyW1Cd');

%==========================================================================
% Cu
%==========================================================================
% Norderney data
NneyW1_Cu = 'NneyW1.xlsx'; 
sheetCu = 2; % sheet containing all information concerning Cu

% possible test if Matlab gets the sheet names
% [status,sheets] = xlsfinfo(NneyW1_Cu )

%NneyW1Cu = xlsread(NneyW1_Cu); % load all data from the excel file

% load the relevant data from the excel file
% dates
x1Range_date = 'H2:H13';
NneyW1Cu_date = xlsread(NneyW1_Cu, sheetCu, x1Range_date);

% delete cells without any data 
NneyW1Cu_date(any(isnan(NneyW1Cu_date),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

% values
x1Range_values = 'P2:P13';
NneyW1Cu_values = xlsread(NneyW1_Cu, sheetCu, x1Range_values);

% Excel and matlab use another system for dates, therfor in ecxel the dates
% have been transformed to a number (use format cells)
% now we need to correct the data in matlab as matlab counts the days since
% the 1 Jänner 0000 and excel counts the days since the 1 Jänner 1900
% therfore it is neccesary to add the 31. Dez 1899 to the imported data
% with te dunction datenum

NneyW1Cu_date = NneyW1Cu_date + datenum('30DEC1899');

%--------------------------------------------------------------------------
% TEST
% test if the date is  transferred correctly from excel to matlab 
% datestr(NneyW1Cu_date(1))
%--------------------------------------------------------------------------

% shows all dates as a string like 21-Feb-2005
NneyW1Cu_dateStr = datestr(NneyW1Cu_date(1:end));

Cu_in_water_unit = 'ug/L';

% create a timeseries 
NneyW1Cu_values(any(isnan(NneyW1Cu_values),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on


%--------------------------------------------------------------------------
% TEST
% plot(NneyW1Cu_date,NneyW1Cu_values)
% length(NneyW1Cd_values) % use this value later for uncertainty analysis
% NneyW1CdH2O = timeseries(NneyW1Cd_values,tspan);
% -------------------------------------------------------------------------

%**************************************************************************
% interpolate data
%**************************************************************************
% interpolate the data so that they can be used in the other programs
% unit in all programs: days
% matlab method: resample(ts, time, interp_method) (interp_method linear by
% default

% TEST
% disp(NneyW1Cd_date)
% disp(NneyW1Cd_values)

NneyW1Cu = interp1(NneyW1Cu_date,NneyW1Cu_values,tnewStress,'pchip'); % before: NneyW1Cupre
% disp('NneyW1Cu')
% disp(length(NneyW1Cu))
% disp(NneyW1Cu)
% NneyW1Cdpr = NneyW1Cupre([100:2096]);
figure(4)
plot(tnewStress,NneyW1Cu,'k',NneyW1Cu_date,NneyW1Cu_values,'ro')
datetick('x','mmm yy','keepticks')
xlabel('month and year')
ylabel('Cu concentration in \mug/L')
title('Cu')

saveas(gcf, fullfile(newFolder, ['CuConc', '.fig']));
    saveas(gcf, fullfile(newFolder, ['CuConcJPG-', '.jpg']));

Cumatfile = 'Cu_in_water.mat'; % create a file for the stressor
save(Cumatfile,'NneyW1Cu_date', 'NneyW1Cu_values', 'Cu_in_water_unit', 'NneyW1Cu') 

% -------------------------------------------------------------------------
%TEST
%shows all dates as a string like 21-Feb-2005
%tnewStr = datestr(tnew);
%disp(tnewStr)
%figure(4)
%plot(tnew,neyW1Cd_values_interp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enable figures if needed
% figure(4)
% % subplot(2,1,1)
% plot(NneyW1Cu_date,NneyW1Cu_values,'o',tnewStress,NneyW1Cu,'r')
% title('data and interpolated data Cu')
% xlabel('time in days')
% ylabel('Cu in ug/L')
% % test 
% % length(NneyW1Cd_values)
% xlabel('time in days')
% ylabel('Cu in ug/L')

% -------------------------------------------

% *************************************************************************
% create a  time series for later use (e.g. manipulation etc.)
%**************************************************************************
NneyW1Cu_values_ts = timeseries(NneyW1Cu_values,NneyW1Cu_date,'name',...
    'NneyW1Cu');

%==========================================================================
% Pb
%==========================================================================
% Norderney data
NneyW1_Pb = 'NneyW1.xlsx'; 
sheetPb = 4; % sheet containing all information concerning Cd

% possible test if Matlab gets the sheet names
% [status,sheets] = xlsfinfo(NneyW1_Cd )

%NneyW1Cd = xlsread(NneyW1_Cd); % load all data from the excel file

% load the relevant data from the excel file
% dates
x1Range_date = 'H2:H13';
NneyW1Pb_date = xlsread(NneyW1_Pb, sheetPb, x1Range_date);

% delete cells without any data 
NneyW1Pb_date(any(isnan(NneyW1Pb_date),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

% values
x1Range_values = 'P2:P13';
NneyW1Pb_values = xlsread(NneyW1_Pb, sheetPb, x1Range_values);

% delete cells without any data 
NneyW1Pb_values(any(isnan(NneyW1Pb_values),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

% Excel and matlab use another system for dates, therfor in ecxel the dates
% have been transformed to a number (use format cells)
% now we need to correct the data in matlab as matlab counts the days since
% the 1 Jänner 0000 and excel counts the days since the 1 Jänner 1900
% therfore it is neccesary to add the 31. Dez 1899 to the imported data
% with te dunction datenum

NneyW1Pb_date = NneyW1Pb_date + datenum('30DEC1899');

%--------------------------------------------------------------------------
% TEST
% test if the date is  transferred correctly from excel to matlab 
% anser should be 21-Feb-2005
% datestr(NneyW1Cd_date(1))
%--------------------------------------------------------------------------

% shows all dates as a string like 21-Feb-2005
NneyW1Pb_dateStr = datestr(NneyW1Pb_date(1:end));


Pb_in_water_unit = 'ug/L';

%--------------------------------------------------------------------------
% TEST
% length(NneyW1Cd_values) % use this value later for uncertainty analysis
% NneyW1CdH2O = timeseries(NneyW1Cd_values,tspan);
% -------------------------------------------------------------------------

%**************************************************************************
% interpolate data
%**************************************************************************
% interpolate the data so that they can be used in the other programs
% unit in all programs: days
% matlab method: resample(ts, time, interp_method) (interp_method linear by
% default

% TEST
% disp(NneyW1Cd_date)
% disp(NneyW1Cd_values)

NneyW1Pb = interp1(NneyW1Pb_date,NneyW1Pb_values,tnewStress,'pchip');

Pbmatfile = 'Pb_in_water.mat'; % create a file for the stressor
save(Pbmatfile,'NneyW1Pb_date', 'NneyW1Pb_values', 'Pb_in_water_unit','NneyW1Pb') 
load(Pbmatfile,'NneyW1Pb_date', 'NneyW1Pb_values', 'Pb_in_water_unit','NneyW1Pb')

figure(5)
plot(tnewStress,NneyW1Pb,'k',NneyW1Pb_date,NneyW1Pb_values,'ro')
datetick('x','mmm yy','keepticks')
xlabel('month and year')
ylabel('Pb concentration in \mug/L')
title('Pb')

saveas(gcf, fullfile(newFolder, ['PbConc', '.fig']));
    saveas(gcf, fullfile(newFolder, ['PbConcJPG-', '.jpg']));

% -------------------------------------------------------------------------
%TEST
%shows all dates as a string like 21-Feb-2005
%tnewStr = datestr(tnew);
%disp(tnewStr)
%figure(5)
%plot(tnew,neyW1Cd_values_interp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enable figures if needed
% enable figures if needed
% figure(5)
% % subplot(2,1,1)
% plot(NneyW1Pb_date,NneyW1Pb_values,'o',tnewStress,NneyW1Pb,'r')
% title('data and interpolated data Pb')
% xlabel('time in days')
% ylabel('Pb in ug/L')
% % test 
% % length(NneyW1Cd_values)
% xlabel('time in days')
% ylabel('Pb in ug/L')

% -------------------------------------------------------------------------

% *************************************************************************
% create a  time series for later use (e.g. manipulation etc.)
%**************************************************************************
NneyW1Pb_values_ts = timeseries(NneyW1Pb_values,NneyW1Pb_date,'name',...
    'NneyW1Pb');

%==========================================================================
% Zn
%==========================================================================

% Norderney data
NneyW1_Zn = 'NneyW1.xlsx'; 
sheetZn = 6; % sheet containing all information concerning Cd

% possible test if Matlab gets the sheet names
% [status,sheets] = xlsfinfo(NneyW1_Cd )

%NneyW1Cd = xlsread(NneyW1_Cd); % load all data from the excel file

% load the relevant data from the excel file
% dates
x1Range_date = 'H2:H13';
NneyW1Zn_date = xlsread(NneyW1_Zn, sheetZn, x1Range_date);

% delete cells without any data 
NneyW1Zn_date(any(isnan(NneyW1Zn_date),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

% delete cells without any data (any(isnan(NneyW1Cd_date),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

% values
x1Range_values = 'P2:P13';
NneyW1Zn_values = xlsread(NneyW1_Zn, sheetZn, x1Range_values);

% delete cells without any data 
NneyW1Zn_values(any(isnan(NneyW1Zn_values),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

% Excel and matlab use another system for dates, therfor in ecxel the dates
% have been transformed to a number (use format cells)
% now we need to correct the data in matlab as matlab counts the days since
% the 1 Jänner 0000 and excel counts the days since the 1 Jänner 1900
% therfore it is neccesary to add the 31. Dez 1899 to the imported data
% with te dunction datenum

NneyW1Zn_date = NneyW1Zn_date + datenum('30DEC1899');

%--------------------------------------------------------------------------
% TEST
% test if the date is  transferred correctly from excel to matlab 
% anser should be 21-Feb-2005
% datestr(NneyW1Cd_date(1))
%--------------------------------------------------------------------------

% shows all dates as a string like 21-Feb-2005
NneyW1Zn_dateStr = datestr(NneyW1Zn_date(1:end));


Zn_in_water_unit = 'ug/L';

%--------------------------------------------------------------------------
% TEST
% length(NneyW1Cd_values) % use this value later for uncertainty analysis
% NneyW1CdH2O = timeseries(NneyW1Cd_values,tspan);
% -------------------------------------------------------------------------

%**************************************************************************
% interpolate data
%**************************************************************************
% interpolate the data so that they can be used in the other programs
% unit in all programs: days
% matlab method: resample(ts, time, interp_method) (interp_method linear by
% default

% TEST
% disp(NneyW1Cd_date)
% disp(NneyW1Cd_values)

NneyW1Zn = interp1(NneyW1Zn_date,NneyW1Zn_values,tnewStress,'pchip');

Znmatfile = 'Zn_in_water.mat'; % create a file for the stressor

save(Znmatfile,'NneyW1Zn_date', 'NneyW1Zn_values', 'Zn_in_water_unit', 'NneyW1Zn','NneyW1Zn')
load(Znmatfile,'NneyW1Zn_date', 'NneyW1Zn_values', 'Zn_in_water_unit', 'NneyW1Zn','NneyW1Zn')

figure(6)
plot(tnewStress,NneyW1Zn,'k',NneyW1Zn_date,NneyW1Zn_values,'ro')
datetick('x','mmm yy','keepticks')
xlabel('month and year')
ylabel('Zn concentration in \mug/L')
title('Zn')

saveas(gcf, fullfile(newFolder, ['ZnConc', '.fig']));
    saveas(gcf, fullfile(newFolder, ['ZnConcJPG-', '.jpg']));

% save NneyW1Zn
% -------------------------------------------------------------------------
%TEST
%shows all dates as a string like 21-Feb-2005
%tnewStr = datestr(tnew);
%disp(tnewStr)
%figure(6)
%plot(tnew,neyW1Cd_values_interp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enable figures if needed
% figure(6)
% % subplot(2,1,1)
% 
% plot(NneyW1Zn_date,NneyW1Zn_values,'o',tnewStress,NneyW1Zn,'r')
% title('data and interpolated data Zn')
% xlabel('time in days')
% ylabel('Zn in ug/L')
% -------------------------------------------------------------------------

% *************************************************************************
% create a  time series for later use (e.g. manipulation etc.)
%**************************************************************************
NneyW1Zn_values_ts = timeseries(NneyW1Zn_values,NneyW1Zn_date,'name',...
    'NneyW1Zn');

%==========================================================================
% Temperature
%==========================================================================

% Norderney data
HAMSOM_Temp = 'HAMSOM_GBtemp.xlsx'; 
sheetTemp = 1; % sheet containing all information concerning Temp

% load the relevant data from the excel file
% dates
x1Range_date = 'D2:D1933';
HAMSOMTemp_date = xlsread(HAMSOM_Temp, sheetTemp, x1Range_date);

% delete cells without any data 
HAMSOMTemp_date(any(isnan(HAMSOMTemp_date),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

% values
x1Range_values = 'E2:E1933';
HAMSOMTemp_values = xlsread(HAMSOM_Temp, sheetTemp, x1Range_values);

% delete cells without any data 
HAMSOMTemp_values(any(isnan(HAMSOMTemp_values),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on


% Excel and matlab use another system for dates, therfor in ecxel the dates
% have been transformed to a number (use format cells)
% now we need to correct the data in matlab as matlab counts the days since
% the 1 Jänner 0000 and excel counts the days since the 1 Jänner 1900
% therfore it is neccesary to add the 31. Dez 1899 to the imported data
% with te dunction datenum

HAMSOMTemp_date = HAMSOMTemp_date + datenum('30DEC1899');

%--------------------------------------------------------------------------
% TEST
% test if the date is  transferred correctly from excel to matlab 
% anser should be 21-Feb-2005
% datestr(NneyW1Cd_date(1))
%--------------------------------------------------------------------------

% shows all dates as a string like 21-Feb-2005
HAMSOMTemp_dateStr = datestr(HAMSOMTemp_date(1:end));


Temp_in_water_unit = '°C';

HAMSOMTemp = HAMSOMTemp_values';
disp('Temp length')
disp(length(HAMSOMTemp))
%--------------------------------------------------------------------------
% TEST
% length(NneyW1Cd_values) % use this value later for uncertainty analysis
% NneyW1CdH2O = timeseries(NneyW1Cd_values,tspan);
% -------------------------------------------------------------------------

%**************************************************************************
% interpolate data
%**************************************************************************
% no interpolation needed as these are already modelled data

% NneyW1Temp = interp1(HAMSOMTemp_date,HAMSOMTemp,tnewStress,'pchip');
Tempmatfile = 'Temp_in_water.mat'; % create a file for the stressor
save(Tempmatfile,'HAMSOMTemp_date', 'HAMSOMTemp', 'Temp_in_water_unit','HAMSOMTemp') 
load(Tempmatfile,'HAMSOMTemp_date', 'HAMSOMTemp', 'Temp_in_water_unit','HAMSOMTemp') 

% save HAMSOMTemp
% -------------------------------------------------------------------------
%TEST
%shows all dates as a string like 21-Feb-2005
%tnewStr = datestr(tnew);
%disp(tnewStr)
%figure(1)
%plot(tnew,neyW1Cd_values_interp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enable figures if needed
% figure(7)
% % subplot(2,1,1)
% plot(tnewStress,HAMSOMTemp,'r')
% title('HAMSOM model data - Temperature')
% xlabel('time in days')

% ylabel('temperature in °C')
figure(7)
plot(HAMSOMTemp_date,HAMSOMTemp_values)
datetick('x','mmm yy','keepticks')
xlabel('month and year')
ylabel('temperature in °C')
title('Temperature')
xlabel('month and year')
ylabel('temperature in °C')

saveas(gcf, fullfile(newFolder, ['Temp', '.fig']));
    saveas(gcf, fullfile(newFolder, ['TempJPG-', '.jpg']));
% -------------------------------------------------------------------------

% *************************************************************************
% create a  time series for later use (e.g. manipulation etc.)
%**************************************************************************
HAMSOMTemp_values_ts = timeseries(HAMSOMTemp_values,HAMSOMTemp_date,...
    'name','HAMSOMTemp');

%==========================================================================
% Salinity
%==========================================================================

% Norderney data
HAMSOM_sal = 'HAMSOM_GBsal.xlsx'; 
sheetsal = 1; % sheet containing all information concerning Cd

% possible test if Matlab gets the sheet names
% [status,sheets] = xlsfinfo(HAMSOM_sal)

% load the relevant data from the excel file
% dates
x1Range_date = 'D2:D1933';
HAMSOMsal_date = xlsread(HAMSOM_sal, sheetsal, x1Range_date);

% delete cells without any data 
HAMSOMsal_date(any(isnan(HAMSOMsal_date),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

% values
x1Range_values = 'E2:E1933';
HAMSOMsal_values = xlsread(HAMSOM_sal, sheetsal, x1Range_values);

% delete cells without any data 
HAMSOMsal_values(any(isnan(HAMSOMsal_values),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on

% Excel and matlab use another system for dates, therfor in ecxel the dates
% have been transformed to a number (use format cells)
% now we need to correct the data in matlab as matlab counts the days since
% the 1 Jänner 0000 and excel counts the days since the 1 Jänner 1900
% therfore it is neccesary to add the 31. Dez 1899 to the imported data
% with te dunction datenum

HAMSOMsal_date = HAMSOMsal_date + datenum('30DEC1899');

%--------------------------------------------------------------------------
% TEST
% test if the date is  transferred correctly from excel to matlab 
% anser should be 21-Feb-2005
% datestr(NneyW1Cd_date(1))
%--------------------------------------------------------------------------

% shows all dates as a string like 21-Feb-2005
HAMSOMsal_dateStr = datestr(HAMSOMsal_date(1:end));


salinity_unit = 'psu';

% create a timeseries 
HAMSOMsal_values(any(isnan(HAMSOMsal_values),2),:) = []; % first remove all 
% values which have not been measured and contain an empty field, which 
% can can create error messages later on
HAMSOMsal = HAMSOMsal_values';
%--------------------------------------------------------------------------
% TEST
% length(NneyW1Cd_values) % use this value later for uncertainty analysis
% NneyW1CdH2O = timeseries(NneyW1Cd_values,tspan);
% -------------------------------------------------------------------------

%**************************************************************************
% interpolate data
%**************************************************************************
% no interpolation needed as these are already modelled data

disp('length salinity')
disp(length(HAMSOMsal))

%NneyW1Sal = interp1(HAMSOMsal_date,HAMSOMsal,tnewStress,'pchip');
Salmatfile = 'salinity.mat'; % create a file for the stressor
save(Salmatfile,'HAMSOMsal_date', 'HAMSOMsal', 'salinity_unit', 'HAMSOMsal') 
load(Salmatfile,'HAMSOMsal_date', 'HAMSOMsal', 'salinity_unit', 'HAMSOMsal') 

save HAMSOMsal
% -------------------------------------------------------------------------
%TEST
%shows all dates as a string like 21-Feb-2005
%tnewStr = datestr(tnew);
%disp(tnewStr)
%figure(1)
%plot(tnew,neyW1Cd_values_interp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% enable figures if needed
% figure(8)
% % subplot(2,1,1)
% plot(tnewStress,HAMSOMsal,'r')
% title('HAMSOM model data - Salinity')
% xlabel('time in days')
% ylabel('salinity')

% plot(HAMSOMsal_date,HAMSOMsal_values)
% title('data')
% xlabel('time in days')
% ylabel('salinity in psu')
% % test 
% % length(NneyW1Cd_values)
% % length(tnew)
% tspan = 1:2096;
% length(tspan)
% subplot(2,1,2)
% plot(tspan,HAMSOMsal_values)
% title('HAMSOM model data')
% xlabel('time in days')
% ylabel('salinity in psu')
% -------------------------------------------------------------------------

% *************************************************************************
% create a  time series for later use (e.g. manipulation etc.)
%**************************************************************************
HAMSOMsal_values_ts = timeseries(HAMSOMsal_values,HAMSOMsal_date,...
    'name','HAMSOMsal');
%==========================================================================
%==========================================================================
EnvData = [NneyW1Cd;NneyW1Cu;NneyW1O2;NneyW1Pb;NneyW1pH;...
    NneyW1Zn;HAMSOMTemp;HAMSOMsal];
end

