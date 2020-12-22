% Script to produce the graphs

function SinglStressorGraphs

opt = load('parametersSingleStressor.mat');

beta = opt.BestParams(1);
mu = opt.BestParams(2);
alfa = opt.BestParams(3);

% plots for showing the model for the eperimental dataset used
BestParams = [beta mu alfa];

N0Exp = 100; % initial poluation size or percent function
TExpx = 0:21; %days
% length(T)

CtabExpy = [0 500 1000 2000 3000 4000 5000 10000 25000]; % concentration eg in ug/L 

for i = 1:length(CtabExpy)
    YExp = SingleStressorModell(CtabExpy(i),BestParams,TExpx,N0Exp);
    NmodExp(:,i) = YExp(:,2);
end

CtabnewXpre = [];
for j = 1:length(CtabExpy)
    CtabnewXAdd = ones(length(TExpx))*CtabExpy(j);
    CtabnewXpre = [CtabnewXpre;CtabnewXAdd];
end

CtabnewX = CtabnewXpre(:,1);

TExpnewpre = [];
for m = 1:length(CtabExpy)
    TExpnewXAdd = TExpx;
    TExpnewpre = [TExpnewpre,TExpnewXAdd];
end

CtabnewX = CtabnewXpre(:,1);
TExpnewY = TExpnewpre';

NewVector = [CtabnewX,TExpnewY];

[NmodExp, emodExp, amodExp] = ModelExp(BestParams,CtabExpy,TExpx,N0Exp);
NmodExpGraph = NmodExp(:);

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
    100 98 97 80 1 1 45 1 9]; % example of a dataset showing responses of the organism for different exposure times and concentrations

NexpGraph = Nexp(:);

figure(1)
[xq,yq] = meshgrid(0:25000/30:25000, 0:21/30:21);
vq = griddata(CtabnewX,TExpnewY,NmodExpGraph,xq,yq);

mesh(xq,yq,vq,'FaceAlpha',0.3)
hold on

% create the graph
plot3(CtabnewX,TExpnewY,NexpGraph,'o')
ax = gca;
ax.FontSize = 15;     % set font size for axes labels
title('model and literature data')
legend('model','sample points','Location','NorthWest')
xlabel('concentration in ug/L','FontSize',15)
ylabel('exposure time in days','FontSize',15)
zlabel('observation','FontSize',15)

YExpz = [];
for j = 1:length(CtabExpy)
    YExpzAdd = Schadstoff_ModellCd_(CtabnewX(j),BestParams,TExpx,N0Exp);
    YExpz = [YExpz;YExpzAdd];
end

function [NmodExp, emodExp, amodExp] = ModelExp(BestParams,CtabExpy,T,N0)
    % N0 is the first value of the matrix to start (position) (e.g. initial population)
    for o = 1:length(CtabExpy)
        YExpOut = Schadstoff_ModellCd_(CtabExpy(o),BestParams,T,N0);
        NmodExp(:,o) = YExpOut(:,2);
        emodExp(:,o) = YExpOut(:,3);
        amodExp(:,o) = YExpOut(:,4);
    end
end


% plots for general model behavior
C = 30; % here in ug/ L,  concentration or intensity of the stressor

parameter0 = [beta,mu,alfa];

xbeta   = [0.2,0.5,1,2]*beta;   
xmu     = [0.1,1,100]*mu;       
xalfa   = [0,5,50,150]*alfa;    
xC      = [0.5,0.2,1,2]*C;     
xpar = {xbeta;xmu;xalfa;xC};
N0   = 100;                     % 100% function, 100% survival or similar

% time steps
T = 0:1:1000;

for par = 1:3                  % loop for the parameters
    figure(par+1); clf; hold on
    parameter = parameter0;
    L = length(xpar{par});
    
    % axis properties
    Cmax = 12;
    mmax = 0.5;
    
    for i = 1:L               	
        if par==1; parameter(1) = xbeta(i);
            Titel = ['effect rate  ','beta =  ',num2str(xbeta)]; end
        if par==2; parameter(2) = xmu(i);
            Titel = ['1/time delay  ','mu =  ',num2str(xmu)]; end
        if par==3; parameter(3) = xalfa(i);
            Titel = ['acclimation  ','alfa =  ',num2str(xalfa)]; end
        if par==7; C = xC(i);
            Titel = ['ext. concentration      ','C =  ',num2str(xC)]; end
        
        % call of the main model
        Y = SSinglStressor_Model(C,parameter,T,N0);
        Cint = Y(:,1);
        N    = Y(:,2);
        a    = Y(:,3);
        m    = Y(:,4);
        
        Cmax = Plot(1,i,T,Cint,'r',Cmax,' ',i==L);
        title(Titel)
        ylabel('internal concentration in ug/g')
        Plot(2,i,T,N,'g',N0,' ',i==L);
        ylabel('percent survival')
        Plot(3,i,T,a,'b',0.02,' ',i==L);
        ylabel('acclimatisation factor')
        %Plot(4,i,T,m,'k',0.07,' ',i==L);
        mmax = Plot(4,i,T,m,'k',mmax,' ',i==L);
        ylabel('effect strength')
        xlabel('time in days')
    end
end
end

% result plot
function Ymax = Plot(sub,i,T,Y,col,Ymax,txt,letzt)
tend = max(T);
if max(Y)>Ymax
    Ymax = max(Y);
end
subplot(4,1,sub);
hold on;
axis([0 tend 0 Ymax])
if letzt
    text(0.65*tend,0.9*Ymax,txt);
    col = [col,'-'];
else
    switch i
        case 1
            % low value
            col = [col,':'];
        case 2
            % moderate value
            col = [col,'-.'];
        case 3
            % high value
            col = [col,'--'];
    end
end
plot(T,Y,col);

end