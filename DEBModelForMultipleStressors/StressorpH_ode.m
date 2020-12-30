% differential equation for the effect of a single stressor
% dY represents dY/dt etc.
function dY = StressorpH_ode(t,Y,pH_)

% Parameter values - from literature

pH0pt = 8.1; %optimum  value based on Heinemann et al. 2012
pHTol = 7.7; % based on Thomsen and Melzner 2010



% paramters calculated (04.06.2018)
betapH = 0.195191063466371;
mupH = 0.001874386292647;
alfapH = 8.544815817237917e-06;


% environmental data
pHW = pH_; 


NpH       = Y(1);
epH       = Y(2);
apH       = Y(3);


%--------------------------------------------------------------------
dN      = -epH*NpH;
e_equpH   = betapH*(abs(pHTol-pHW))/(abs(pH0pt-pHTol))*(1-apH);
depH      = mupH*(e_equpH-epH);

dapH      = (alfapH*epH)*(1-apH);

dY      = [dN;depH;dapH];

%----------------------------------------------------------------------

