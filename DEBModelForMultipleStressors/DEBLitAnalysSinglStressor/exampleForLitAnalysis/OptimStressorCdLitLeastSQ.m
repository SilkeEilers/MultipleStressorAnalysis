% stressor model based on Kooijman and Bedaux 1996 Analysis of toxicity
% data page 64 "Toxicity test on survival". Model modified. 
% Data for the blue mussel Mytilus edulis

saveparametersCd = 'parametersCd.mat';

% Literature data
% mortality data derived from graphs by Sunila 1981 Toxicity of copper and
% cadmium to Mytilus edulis in brakish water Ann. Zoo. Fennici 18. 213-223


Nexp = [100 100 100 100 100 100 100 100 100; 
    100 100 100 100 98 92 97 92 67; 
    100 100 100 100 93 91 96 91 65; 
    100 100 100 100 91 90 96 90 64; 
    100 100 99 100 88 85 96 86 63;
    100 100 100 100 100 98 98 98 88; 
    100 100 100 100 100 94 97 93 71; 
    100 100 99 99 85 81 95 83 62; 
    100 100 99 99 82 80 95 81 62; 
    100 100 99 98 80 78 95 81 61; 
    100 100 99 97 60 67 95 78 61; 
    100 100 99 97 47 55 93 67 60;
    100 100 99 97 21 31 92 50 60;
    100 100 99 97 12 24 92 41 58;
    100 100 99 90 10 11 90 27 52;
    100 100 98 85 9 5 86 18 47;
    100 100 98 85 7 3 85 11 40;
    100 99 99 82 5 3 80 8 28;
    100 98 99 82 4 2 68 6 13;
    100 98 98 81 2 2 65 2 11;
    100 98 98 81 2 2 53 2 10;
    100 98 97 80 1 1 45 1 9];


N0 = Nexp(1,1);                 

% data from the experimental setup (conditions)
% exposure times
T = 0:21; %days, data by Sunila et al 1981

% concentrations tested
Ctab = [0 500 1000 2000 3000 4000 5000 10000 25000]; % ug/L calculated from
% ppm unit data by Sunila et al 1981, Fig. 2

% parameters - initial guesses
betaCdIni = 0.002;       % mortality rate pro CCdT, 
                       % dimension 1/(time*conc) per time and int. conc.
muCdIni = 0.002;           % 1/latency time, can have an initity value, then m=mort
alfaCdIni = 0.002;      % speed of adaptation, relative to mortality
ParamsIni = [betaCdIni,muCdIni,alfaCdIni];

% Minimize the difference between the model and the data
%BestParams = lsqnonlin(@(param) FitModel(param,Ctab,T,N0,Nexp),ParamsIni,[0 0 0],[1 50 1]);

%BestParams = fmincon(@(param) FitModel(param,Ctab,T,N0,Nexp),ParamsIni,[],[],[],[],[0 0 0],[1 50 1]);
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(param) FitModel(param,Ctab,T,N0,Nexp),ParamsIni,[],[],[],[],[0 0 0],[1 50 1]);
% plot the original and experimental data
%%%ModelNew = Schadstoff_ModellCd_(Ctab(i),BestParams,T,N0);


% calculate the standard error
BestParams=x;

% calcualte the model with the best parameters
NmodDecimal = Model(BestParams,Ctab,T,N0);
DiffModelDataDecimal = Nexp(:)-NmodDecimal(:);
relAbsDiff = abs(DiffModelDataDecimal./(Nexp(:)+NmodDecimal(:)));

n = length(Nexp(:));

% Standard Error des besten Modells
MeanError = mean(abs(DiffModelDataDecimal));
MeanRelError = mean(relAbsDiff);
StandardDeviation = std(abs(DiffModelDataDecimal));
StandardRelDeviation = std(relAbsDiff);
StandardError = std(abs(DiffModelDataDecimal))/sqrt(length(DiffModelDataDecimal));

% old (only for comparison)
% initError = sum(FitModel(ParamsIni,Ctab,T,N0,Nexp).^2)
% bestError = sum(FitModel(BestParams,Ctab,T,N0,Nexp).^2)


% display the results
disp('Nexp')
disp(Nexp)
disp('Nmod')
disp(round(Model(BestParams,Ctab,T,N0)))
disp('Diff')
disp(Nexp-round(Model(BestParams,Ctab,T,N0)))

% optimized parameters
disp('best params')
disp(BestParams)

% mean error
disp('mean error')
disp(MeanError)
disp('mean relative error')
disp(MeanRelError)

% standard deviation
disp('standard deviation')
disp(StandardDeviation)
disp('standard relative deviation')
disp(StandardRelDeviation)

% standard error
disp('standard error')
disp(StandardError)


% number of iterations
disp('number of iterations')
disp(output.iterations)

% save parameters
save(saveparametersCd,'BestParams')

% functions
function DiffModelData = FitModel(Params,Ctab,T,N0,Nexp)
    Nmod=Model(Params,Ctab,T,N0);
    DiffModelDatapre = Nexp(:)-Nmod(:);
    DiffModelData=sum(abs(DiffModelDatapre));
    %DiffModelData=sum((DiffModelData).^2);
end

function Nmod = Model(Params,Ctab,T,N0)
    % N0 is the first value of the matrix to start (position) (e.g. initial population)
    for i = 1:length(Ctab)
        Y = Stressor_modelCd_(Ctab(i),Params,T,N0);
        Nmod(:,i) = Y(:,2);
    end
end



