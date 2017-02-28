function [wat_cont,new_WTD,new_peat_wat,ET,ET_temp_fact,transmis,Roff,Ron,flag] = hpm_WatBal12pf(WTD,DEL_WAT,PEATWATER,DENS,THICK,DEPTH,POROSITY,TOT_porosity,Zstar,Zbottom,ONEVEC,ZEROVEC,params,time,ALT,Tmean)

% v.12pf -- no change from v.10pf

% v.10pf adds constraint of active layer depth for transmissivity

% v.9 includes change in runoff to start peat in low areas

% no change from v.7

% function calculates new WTD and peat water content, ET, Transmissivity,Run-on and Run-off

% v.7 finds the new water table for net addition or loss of water

% OUTPUT VARIABLES
% ------------------
% wat_cont is cohort WFPS*porosity above the water table (m3/m3); below WT wat_cont = porosity
% new_WTD is final water table depth below surface [m]
% new_peat_water is final water content of peat column (saturated + unsaturated) [m or m3/m2]
% ET is annualized ET [m/y], a function of WTD
% transmis is relative tranmissivity [-], a funciton of hydraulic conductivity profile and WTD; hydr. cond. profile a function of peat density (degree of decomposition)
% Roff is annualized runoff [m/y], a function of WTD and transmissivity
% Ron is annualized runon [m/y] a function of WTD and peat depth????

% INPUT VARIABLES
% ------------------
% WTD is initial water table depth below surface [m]
% DEL_WAT is net water added to peat profile in timestep [m or m3/m2]
% PEATWATER is initial total water content of peat column (saturated + unsaturated) [m or m3/m2]
% DENS is cohort bulk density profile [kg/m3]
% THICK is vector of cohort thicknesses (m)
% DEPTH is vector of cohort mid-point depths (m, positive down)
% POROSITY is vector of cohort porosities (- or m3/m3)
% TOT_porosity is total peat column pore volume/area (m or m3/m2)
% Tmean is mean annual temperature (for ET factor)
% Zstar is factor for unsaturated WFPS, and is function of bulk density
% Zbottom is array of depths (m) from top of peat to bottom of cohort
% ONEVEC & SEROVEC are vectors of ones and zeros
% params is model input parameters
% time is previous year's model time step (years) - used as a counter
% ALT is current year's active layer depth (as water balance is calculated before active layer depth)
% Tmean is mean annual temperature
% ------------------

zwtd = zeros(params.sim_len,1);  % trying to fix a crash that doesn't make sense; assumes annual time step
flag = 0;

%% ******************
% *** LOCATE WATER TABLE ***
% ******************

total_water = PEATWATER + max(0,-WTD);  % total water in peat profile [m]

new_water = total_water + DEL_WAT;


if (max(Zbottom) > ALT)   % only check for WTD below ALT if peat is thicker than ALT
    counter_pf = find(Zbottom > ALT, 1);
else
    counter_pf = time-5;   % move it up a bit for 'safety'
end

zwtd_pf = DEPTH - sum(THICK(1:counter_pf));
zwtd_pf = max(ZEROVEC, -zwtd_pf);
wat_cont_pf = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd_pf ./ Zstar);  
peat_water_pf = sum((wat_cont_pf .* THICK) .* POROSITY);

if (new_water >= TOT_porosity)  % peat is saturated, WTD at or above surface (i.e., WTD <= 0)
    
    new_WTD = -(new_water - TOT_porosity);
    wat_cont = ONEVEC;
    new_peat_wat = TOT_porosity;
    counter = 1;
    
elseif (WTD > max(Zbottom))  % put WT at bottom of peat and calculate final layer water contents and total peat water for new WT
    
    new_WTD = max(Zbottom) - 0.20 * params.start_depth  % move WTD 20% of 'start_depth' above base of peat profile 
    peat_ht = max(Zbottom)
    flag = 1; 
    zwtd = DEPTH - new_WTD;
    zwtd = max(ZEROVEC, -zwtd);
    wat_cont = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd./Zstar);
    new_peat_wat = sum((wat_cont .* THICK) .* POROSITY);
    counter = time;
    
elseif (new_water <= peat_water_pf)  % prevent WT from going deeper than active layer
    new_WTD = ALT;
    zwtd = DEPTH - new_WTD;
    zwtd = max(ZEROVEC, -zwtd);
    wat_cont = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd./Zstar);
    new_peat_wat = sum((wat_cont .* THICK) .* POROSITY);
    counter = counter_pf;

else

    counter = find(Zbottom > WTD, 1);  % layer # with previous water table
    
    if (counter > 1.5)
        zwtd_tmp = (DEPTH) - sum(THICK(1:counter));
        zwtd_tmp = max(ZEROVEC, -zwtd_tmp);
        wat_cont_tmp = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd_tmp./Zstar);   % see notes and file 'anoxia & bulk dens & WFPS % profile.xls')
        peat_water_bot = sum((wat_cont_tmp .* THICK) .* POROSITY);
        zwtd_tmp = (DEPTH) - sum(THICK(1:counter-1));
        zwtd_tmp = max(ZEROVEC, -zwtd_tmp);
        wat_cont_tmp = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd_tmp./Zstar);   % see notes and file 'anoxia & bulk dens & WFPS % profile.xls')
        peat_water_top = sum((wat_cont_tmp .* THICK) .* POROSITY);
    else
        zwtd_tmp = (DEPTH) - sum(THICK(1:counter));
        zwtd_tmp = max(ZEROVEC, -zwtd_tmp);
        wat_cont_tmp = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd_tmp./Zstar);   % see notes and file 'anoxia & bulk dens & WFPS % profile.xls')
        peat_water_bot = sum((wat_cont_tmp .* THICK) .* POROSITY);
        peat_water_top = TOT_porosity;
    end
    
    if (DEL_WAT < 0)    % WTD is increasing (water table is dropping)
    
        while (peat_water_bot > new_water)   % keep moving deeper, a layer at a time, until peat water < new water
            counter = counter + 1;
            if (counter > time + 0.5)   % don't let WT go below peat into sub-soil
                counter = counter - 1;
%                zwtd_tmp = (DEPTH - THICK/2) - sum(THICK(1:counter));
                zwtd_tmp = DEPTH - sum(THICK(1:counter));
                zwtd_tmp = max(ZEROVEC, -zwtd_tmp);
                wat_cont_tmp = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd_tmp./Zstar); 
                peat_water_bot = sum((wat_cont_tmp .* THICK) .* POROSITY);
%                new_water = peat_water_tmp;            
                counter1 = counter - 1;
%                zwtd_tmp = (DEPTH - THICK/2) - sum(THICK(1:counter));
                zwtd_tmp = DEPTH - sum(THICK(1:counter1));
                zwtd_tmp = max(ZEROVEC, -zwtd_tmp);
                wat_cont_tmp = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd_tmp./Zstar); 
                peat_water_top = sum((wat_cont_tmp .* THICK) .* POROSITY);
                break;

            else
                peat_water_top = peat_water_bot;
            end
            zwtd_tmp = DEPTH - sum(THICK(1:counter));
            zwtd_tmp = max(ZEROVEC, -zwtd_tmp);
            wat_cont_tmp = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd_tmp./Zstar); 
            peat_water_bot = sum((wat_cont_tmp .* THICK) .* POROSITY);
        end
        
    else         % DEL_WAT > 0 so WTD is decreasing (water table is rising)
        
        while (peat_water_top < new_water)   % keep moving shallower, a layer at a time, until peat water > new water
            counter = counter - 1;
            if (counter < 1.5)  % WT must be in top layer
                counter = 1;
%                zwtd_tmp = (DEPTH - THICK/2) - sum(THICK(1:counter));
                zwtd_tmp = DEPTH - THICK(1);
                zwtd_tmp = max(ZEROVEC, -zwtd_tmp);
                wat_cont_tmp = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd_tmp./Zstar); 
                peat_water_bot = sum((wat_cont_tmp .* THICK) .* POROSITY);
                peat_water_top = TOT_porosity;
                break;

            else
                peat_water_bot = peat_water_top;
            end
            
%            zwtd_tmp = (DEPTH - THICK/2) - sum(THICK(1:counter));
            zwtd_tmp = DEPTH - sum(THICK(1:counter-1));
            zwtd_tmp = max(ZEROVEC, -zwtd_tmp);
            wat_cont = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd_tmp./Zstar);
            peat_water_top  = sum((wat_cont .* THICK) .* POROSITY);

        end
        
    end
    
    new_WTD = sum(THICK(1:counter)) - THICK(counter) * (new_water - peat_water_bot)/(eps + peat_water_top - peat_water_bot);

% calculate final layer water contents and total peat water for new WT
    zwtd = DEPTH - new_WTD;
    zwtd = max(ZEROVEC, -zwtd);
    wat_cont = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd./Zstar);
    new_peat_wat = sum((wat_cont .* THICK) .* POROSITY);
        
end
        
%% ******************
% *** EVAPOTRANSPIRATION ***
% ******************

% function calculates annual ET (m)
% follows revised version of PAM from Roulet
% reference for PAM is:
% Hilbert D, NT Roulet, TR Moore. 2000. Modelling and analysis of peatlands as dynamical systems, J. Ecology. 88:230-242.

if (WTD < params.ET_wtd_1)
    ET = params.ET_0;
elseif (WTD <= params.ET_wtd_2)
    ET = params.ET_0 / (1 + (params.ET_wtd_3-1)*params.ET_param * (WTD - params.ET_wtd_1));
else
    ET = params.ET_0 / params.ET_wtd_3;
end

% adjust ET for annual temperature
%    5% change in ET per degree C (Brummer et al. 2011; Ag. For. Met. 153, 14-30)

ET_temp_fact = 1 + params.ET_temp_sens * (Tmean - params.ann_temp);
ET = ET * ET_temp_fact;

%% ******************
% *** RUN-OFF ***
% ******************

% function calculates peat hydraulic transmissivity 
%    first calculates hydraulic conductivity of cohorts down profile
%    hydraulic conductivity as function of bulk density from Radforth (1977)

% Calculate hydraulic conductivity down the cohort profile
% log_10(K) = (150 - 3*rho)/70  From Radforth (1977)
% 150/70 = 2.142857 ; 3/70 = 0.042857; 10^(150/70) = 138.94955 ; 10^(-3/70) = 0.906031

trans_counter1 = find(DEPTH > new_WTD,1);

trans_counter_pf = find(DEPTH > ALT,1);

trans_counter = min(trans_counter1, trans_counter_pf);

% counter2 = find(DEPTH > parameter,1);  % parameter = depth above which water table rarely goes

Log10HydrCond = 2.14287 * ONEVEC - 0.042857 * DENS;
HydrCond = exp(log(10)*Log10HydrCond);

if (trans_counter_pf >= trans_counter1)
%    transmis = sum(THICK(trans_counter1:trans_counter_pf) .* HydrCond(trans_counter1:trans_counter_pf)) / sum(THICK .* HydrCond);    
    transmis = sum(THICK(trans_counter1:trans_counter_pf) .* HydrCond(trans_counter1:trans_counter_pf)) / ...
                   sum(THICK(1:trans_counter_pf) .* HydrCond(1:trans_counter_pf));    
else
    transmis = 0;
end
    
% transmis = sum(THICK(trans_counter:end) .* HydrCond(trans_counter:end)) / sum(THICK(counter2:end) .* HydrCond(counter2:end));

transmis = params.Roff_c3 + (1 - params.Roff_c3) * transmis;  % scale transmissivity range from 0.5 to 1.0

% transmis = 0.92;

% HydrCond = 138.94955 * (ONEVEC * 0.906031) .^ DENSITY;
% HydrCond = 0.01 * ONEVEC;

% calculate relative hydraulic transmissivity

% ALTERNATIVE METHOD #2
% below = DEPTH > WTD;  % creates 'logical' vector with ones for elements with depth>WTD, zeros otherwise
% transmissivity = sum(THICK .* (HydrCond .* below)) / sum(THICK .* HydrCond);

% ALTERNATIVE METHOD #1
% NOTE: this first one gives different results--why?
% THICK1 = THICK;
% THICK1(THICK1<=WTD) = 0;
% transmissivity = sum(THICK1 .* HydrCond) / sum(THICK .* HydrCond);

% function calculates annual water runoff (m/y)
% follows revised version of PAM from Roulet (see ref above)
% reference for PAM is:

runoff1 = params.Roff_c1 * (1 + params.Roff_c2 * (sum(THICK) - params.Roff_c2a));  % modified run-off (March 2010)

if (WTD <= params.Roff_c4)
%    Roff = params.Roff_c4 - WTD + runoff1;
%    Roff = transmis * runoff1 * (1 + (WTD - params.Roff_c4) / (-0.1 + 2 * params.Roff_c4));
   Roff = transmis * runoff1 * min(4,abs(1 + (WTD - params.Roff_c4) / (-0.1 + 2 * params.Roff_c4))); % This removing excess water from site
%    resutling in flooding and NPP crashing.
    Roff  = - WTD;
else
    Roff = transmis * runoff1;
end

%Roff = Roff * 0.75;

Roff = max(transmis * runoff1, -WTD);

% add all water above 0.05 m inundation (WTD = -0.05 m) to runoff
% NOW DOING THIS IN MAIN CODE AFTER STEPPING THROUGH FRACTIONS OF YEAR

%if (new_WTD + Roff < -0.05)
%    Roff = Roff + (-0.05) - (new_WTD + Roff)
%    new_WTD = -0.05;
%end
    
%% ******************
% *** RUN-ON ***
% ******************

% function calculates annual water runon (m/y)
% runon a function of total peat height (just a placeholder for now)

peatdepth = sum(THICK);

if (params.runon_c3 > 0)
    Ron = (params.runon_c3) * (1 - 0.5 * ( 1 + erf((peatdepth - params.runon_c1)/(sqrt(2)*params.runon_c2))));
end

% TESTING A NEW IDEA (APRIL 2009)
% run-on declines with decreasing WTD (i.e., as wter table rises) 

if (params.runon_c3 > 0)
    runon_c4 = 0.0; % WTD depth above which run-on is zero
    runon_c5 = 0.1; % WTD depth below which run-on is maximum
    Ron = Ron * min(1.0, max(0.0, (WTD-runon_c4)/runon_c5));  % Ron factor rises from zero at WTD = 0 to one at WTD = 0.1 m
else
    Ron = 0;
end
return