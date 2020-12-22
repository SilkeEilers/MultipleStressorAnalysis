% differential equation for the effect of a single stressor
% dY represents dY/dt etc.
function dY = StressorTemp_ode(t,Y,Tact)

% Parameter values - from literature

Temp0pt = 15.8; %optimum  value based on Lauzon-Guay et al 2006


TempTolmax = 25; % based on Zittier et al 2015 (MO2, oxygen coinsumption data) 
TempTolmin = -15; % no literature value was found for the lower limit as
                  % a temperature value. However it is known that several
                  % years of ice winter after each other cause effects on
                  % the mussel populations (Witte et al. 2013). Therefore 
                  % a value was chosen, which is typical for ice winters
                  % but this is only an approximation chosen by a best
                  % guess!



% paramters calculated (04.06.2018)
betaTemp = 0.044340766222543;
muTemp = 45.680393049479860; 
alfaTemp = 0.996632117258138; 

% environmental data
Tempglob = Tact; 

NTemp       = Y(1);
eTemp       = Y(2);
aTemp       = Y(3);


%--------------------------------------------------------------------

dN      = - eTemp*NTemp; 
% the difference between Tempglob and TempTol needs to be higher than 0 to
% cause an effect, otherwiae e_equTemp is 0
if Tempglob < 15
    TempTol = TempTolmin;
    %disp('TempTolmin')
else 
    TempTol = TempTolmax;
    %disp('TempTolmax')
end
    
e_equTemp   = betaTemp*(abs(Tempglob-Temp0pt))/(abs(TempTol-Temp0pt))*(1-aTemp);

deTemp      = muTemp*(e_equTemp-eTemp);

daTemp      = alfaTemp*eTemp*(1-aTemp);


dY      = [dN;deTemp;daTemp];
%----------------------------------------------------------------------

