% start value analysis
saveStartValueParamsCd = 'StartValueParamsCd.mat';

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

BestParamValues = [];

% test different start values
% StartValuesB = [0.1 0.2 0.3]; %linspace(0,1);
% StartValuesMu = [0.1 0.2 0.3]; %linspace(0,50);
% StartValuesA = [0.1 0.2 0.3];%linspace(0,1);

StartValuesB = linspace(0,1,10); % numbers between 0 an 1 with 10 points
StartValuesMu = linspace(0,50,10); % numbers between 0 an 50 with 10 points
StartValuesA = linspace(0,1,10); % numbers between 0 an 1 with 10 points

[cb, cmu, ca ] = ndgrid(StartValuesB, StartValuesMu, StartValuesA);
combsStartValues = [cb(:), cmu(:), ca(:)];

% iterate in a loop for different possible combinations of start values

for i = 1:length(combsStartValues)
    
    ParamsIni = combsStartValues(i,:);
    
    [x,fval,exitflag,output,lambda,grad,hessian] = ...
        fmincon(@(param) FitModel(param,Ctab,T,N0,Nexp),...
        ParamsIni,[],[],[],[],[0 0 0],[1 50 1]);
    BestParamValuesAdd = x;
    BestParamValues = [BestParamValues;BestParamValuesAdd];
    disp(i)
end
         

% extract the results for each parameter
BestParamValuesB = BestParamValues(:,1);
BestParamValuesMu = BestParamValues(:,2);
BestParamValuesA = BestParamValues(:,3);

% plot the results
figure(1)
histogram(BestParamValuesB)
title('Best parameter values of beta')

figure(2)
histogram(BestParamValuesMu)
title('Best parameter values of mu')

figure(3)
histogram(BestParamValuesA)
title('Best parameter values of alfa')

% estimate the the most frequent of each parameter with KernelSmoothing
% for each parameter

% beta
[N,edges] = histcounts(BestParamValuesB, 'Normalization', 'probability');
maxProbability = max(N);
PositionMaxProb = find(max(N));
MaxProbValueBetaMinBorder = edges(PositionMaxProb);
MaxProbValueBetaMaxBorder = edges(PositionMaxProb+1);

StartValueBeta = (MaxProbValueBetaMinBorder+MaxProbValueBetaMaxBorder)/2;
% result 0.025

% mu
[N,edges] = histcounts(BestParamValuesMu, 'Normalization', 'probability');
maxProbability = max(N);
PositionMaxProb = find(max(N));
MaxProbValueMuMinBorder = edges(PositionMaxProb);
MaxProbValueMuMaxBorder = edges(PositionMaxProb+1);

StartValueMu = (MaxProbValueMuMinBorder+MaxProbValueMuMaxBorder)/2;
% result: 2.5

% alfa
[N,edges] = histcounts(BestParamValuesA, 'Normalization', 'probability');
maxProbability = max(N);
PositionMaxProb = find(max(N));
MaxProbValueAMinBorder = edges(PositionMaxProb);
MaxProbValueAMaxBorder = edges(PositionMaxProb+1);

StartValueA = (MaxProbValueAMinBorder+MaxProbValueAMaxBorder)/2;
% result: 0.0025

FinalStartValues = [StartValueBeta StartValueMu StartValueA];

save(saveStartValueParamsCd,'BestParamValues', 'FinalStartValues')


toc
t = toc;
save(saveStartValueParametersCd,'BestParamValues')

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
        Y = Schadstoff_ModellCd_(Ctab(i),Params,T,N0);
        Nmod(:,i) = Y(:,2);
    end
end



