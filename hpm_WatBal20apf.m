function [wat_cont,new_WTD,new_PEAT_wat] = hpm_WatBal12a(WTD,del_PeatHt,PEATWATER,TOT_water,DENS,THICK,DEPTH,POROSITY,TOT_porosity,Zstar,Zbottom, ONEVEC,ZEROVEC,params,time,ALT)

% no change from v.8
% no change from v.7a

% function similar to hpm_WatBal6

% v.7a caculates the change in water table position after peat accumulation/loss 
% it doesn't calculate ET, Run-on, Run-off or have and delta water stuff.

% OUTPUT VARIABLES
% ------------------
% wat_cont is cohort WFPS*porosity above the water table (m3/m3); below WT wat_cont = porosity
% new_WTD is final water table depth below surface [m]
% new_peat_water is final water content of peat column (saturated + unsaturated) [m or m3/m2]

% INPUT VARIABLES
% ------------------
% WTD is initial water table depth below surface [m]
% del_PeatHt is net change in total peat column height during year [m]
% PEATWATER is initial total water content of peat column (saturated + unsaturated) [m or m3/m2]
% DENS is cohort bulk density profile [kg/m3]
% THICK is vector of cohort thicknesses (m)
% DEPTH is vector of cohort mid-point depths (m, positive down)
% POROSITY is vector of cohort porosities (- or m3/m3)
% TOT_porosity is total peat column porosity  (- or m3/m3)
% Zstar is factor for unsaturated WFPS, and is function of bulk density
% Zbottom is array of depths (m) from top of peat to bottom of cohort 
% ONEVEC & SEROVEC are vectors of ones and zeros
% params is model input parameters
% time is current model time step (years)
% ------------------

if (max(Zbottom) > ALT)   % only check for WTD below ALT if peat is thicker than ALT
    counter_pf = find(Zbottom > ALT, 1);
else
    counter_pf = time-5;   % move it up a bit for 'safety'
end

zwtd_pf = DEPTH - sum(THICK(1:counter_pf));
zwtd_pf = max(ZEROVEC, -zwtd_pf);
wat_cont_pf = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd_pf ./ Zstar);  
peat_water_pf = sum((wat_cont_pf .* THICK) .* POROSITY);

if (TOT_water >= TOT_porosity)  % peat is saturated, WTD at or above surface (i.e., WTD <= 0)
    
    new_WTD = -(TOT_water - TOT_porosity);
    wat_cont = ONEVEC;
    new_PEAT_wat = TOT_porosity;
    
elseif (WTD > max(Zbottom))  % put WT at bottom of peat and calculate final layer water contents and total peat water for new WT
    
    new_WTD = max(Zbottom);
    zwtd = DEPTH - new_WTD;
    zwtd = max(ZEROVEC, -zwtd);
    wat_cont = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd./Zstar);
    new_PEAT_wat = sum((wat_cont .* THICK) .* POROSITY);

elseif (PEATWATER <= peat_water_pf)  % prevent WT from going deeper than active layer
    new_WTD = ALT;
    zwtd = DEPTH - new_WTD;
    zwtd = max(ZEROVEC, -zwtd);
    wat_cont = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd./Zstar);
    new_PEAT_wat = sum((wat_cont .* THICK) .* POROSITY);

else

    counter = find(Zbottom > WTD, 1);  % current layer # with previous water table
    
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
    
    
    if (peat_water_bot > TOT_water)  % WT is in cohort lower than 'counter'
    
        while (peat_water_bot > TOT_water)   % keep moving deeper, a layer at a time, until peat water < new water
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
        

    else  % peat_water_bot < TOT_Water, so WT is in 'counter' cohort or closer to surface.
    
        while (peat_water_top < TOT_water)   % keep moving up, a layer at a time, until peat water < TOT water
            counter = counter - 1;
            if (counter < 1.5)      % must be in top layer (as it cannot be above peat surface
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
    
    new_WTD = sum(THICK(1:counter)) - THICK(counter) * (TOT_water - peat_water_bot)/(eps + peat_water_top - peat_water_bot);

% calculate final layer water contents and total peat water for new WT
    zwtd = DEPTH - new_WTD;
    zwtd = max(ZEROVEC, -zwtd);
    wat_cont = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd./Zstar);
    new_PEAT_wat = sum((wat_cont .* THICK) .* POROSITY);
        
end

return