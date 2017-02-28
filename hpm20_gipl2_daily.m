function [soil_layer_temp_month soil_node_temp_month swe_day_prev soilTemp ALD1 ALD2 ALD3 max_snow_depth] = ...
    hpm20_gipl2_daily(iyear, T_init_gipl, air_temp_month, monthly_precip, swe_day_prev, depth, thick, annWFPS, poros, dens, params, params_gipl)
    
% uses UAF-GI GIPL2 model to compute daily soil temperature profiles, then averages to monthly
% GIPL2 code supplied by Sergei Marchenko (UAF), worked with version Claire Treat used in her MS
% 

% monthly_precip(1:12) = avg_monthly_precip * (1 + 0.4 * (0.5 - rand(1,12)));
%  mid_month_doy = int32(15 + (30.5 * [0:11]));
mid_month_doy = 15 + (30.5 * [0:11]);

% *****OLD STUFF *********************
% monthly_air_temp_time_series_linear = zeros(num_years *12,4);
% monthly_precip_time_series_linear = zeros(num_years *12,4);
% 
% for i = 1:1:num_years
%     for j = 1:1:12
%         imonth = (i-1)*12 + j;
%         monthly_air_temp_time_series_linear(imonth,1) = monthly_weather_time_series.data(i,1);
%         monthly_air_temp_time_series_linear(imonth,2) = imonth;
%         monthly_air_temp_time_series_linear(imonth,3) = int32(15 + (30.5 * imonth));
%         monthly_air_temp_time_series_linear(imonth,4) = monthly_weather_time_series.data(i,j+1);
%         
%         monthly_precip_time_series_linear(imonth,1) = monthly_weather_time_series.data(i,1);
%         monthly_precip_time_series_linear(imonth,2) = imonth;
%         monthly_precip_time_series_linear(imonth,3) = int32(15 + (30.5 * imonth));
%         monthly_precip_time_series_linear(imonth,4) = monthly_weather_time_series.data(i,j+13)/1000;  % convert mm/month to m/month
%     end
% end

%  --- convert monthly to daily by smooth interpolation
%   Note: so far, precipitation is constant (ann_precip / 366 per day) -- this is only used here for
%   the snowpack

num_months = 12;
num_days = 12 * 30.5;
day_vec = 1:1:num_days;
% day_vec = day_vec';

daily_air_temp_time_series_linear = interp1(mid_month_doy,air_temp_month',day_vec,'pchip','extrap');

%  REMOVE THIS  REMOVE THIS  REMOVE THIS  REMOVE THIS  REMOVE THIS  REMOVE THIS  REMOVE THIS  REMOVE THIS 
%  --- DEBUGGING
%     daily_air_temp_time_series_linear = daily_air_temp_time_series_linear - 2.;
%  REMOVE THIS  REMOVE THIS  REMOVE THIS  REMOVE THIS  REMOVE THIS  REMOVE THIS  REMOVE THIS  REMOVE THIS 

daily_precip_time_series_linear = interp1(mid_month_doy,monthly_precip,day_vec,'pchip','extrap')/30.5;

% placeholder peat values to get code running:
% depth = 0.01 : 0.02 : 2;  % 2 m of peat
totalPeatDepth = max(depth);
% depth_len = length(depth);
% annWFPS = 0.65 * ones(1,depth_len);
% annWFPS(25:end) = annWFPS(25:end) + 0.35;
% annWFPS(1:10) = annWFPS(1:10) - 0.35;
% poros = 0.95 * ones(1,depth_len);

daily_soil_temp_out = zeros(params_gipl.NumberOfSoilComputationNodes,num_days);
monthly_soil_temps = zeros(params_gipl.NumberOfSoilComputationNodes,12);
monthly_layer_temps = zeros(params_gipl.ndepth-1,12);
daily_soil_temps = zeros(params_gipl.NumberOfSoilComputationNodes,366);

% soil_temp_out = zeros(NumberOfSoilComputationNodes,num_months);

hpm_depth_index = zeros(1,params_gipl.NumberOfSoilComputationNodes);
fpeat = zeros(1,params_gipl.NumberOfSoilComputationNodes);
fwat = zeros(1,params_gipl.NumberOfSoilComputationNodes);
fmin = zeros(1,params_gipl.NumberOfSoilComputationNodes);
fporos = zeros(1,params_gipl.NumberOfSoilComputationNodes);
fair = zeros(1,params_gipl.NumberOfSoilComputationNodes);
mean_ann_wfps = zeros(1,params_gipl.NumberOfSoilComputationNodes);
Cond_Th = zeros(1,params_gipl.NumberOfSoilComputationNodes);
Cond_Fr = zeros(1,params_gipl.NumberOfSoilComputationNodes);
Cvol_Th = zeros(1,params_gipl.NumberOfSoilComputationNodes);
Cvol_Fr = zeros(1,params_gipl.NumberOfSoilComputationNodes);
Cvol_Sl = zeros(1,params_gipl.NumberOfSoilComputationNodes);

U1 = zeros(1,params_gipl.NumberOfSoilComputationNodes);
P1 = zeros(1,params_gipl.NumberOfSoilComputationNodes);
Q1 = zeros(1,params_gipl.NumberOfSoilComputationNodes);

% T_air = zeros(1,12);
% P = zeros(1,12);
% snowfall = zeros(1,12);
% rainfall = zeros(1,12);
% snowsublimation = zeros(1,12);
% snowmelt = zeros(1,12);
% swe = zeros(1,12);
% snowDepth = zeros(1,12);
% ALFA = zeros(1,12);

T_air = zeros(1,366);
P = zeros(1,366);
snowfall = zeros(1,366);
rainfall = zeros(1,366);
snowsublimation = zeros(1,366);
snowmelt = zeros(1,366);
swe = zeros(1,366);
snowDepth = zeros(1,366);
ALFA = zeros(1,366);
zerovec = zeros(1,366);

% ALD1 = zeros(1,num_years);  % max annual depth of fully thawed layers 0°C)
% ALD2 = zeros(1,num_years);  % max annual depth of freeze/thaw mid-point (Tfr param)
% ALD3 = zeros(1,num_years);  % min annual depth of never thawed layers (Tfr - FIT params)

    % see Dominik's paper for the unfrozen water function (negligible for peat, so maybe not a huge issue)

% **************************************
%  Inititialize peat and soil thermal properties each year
% **************************************

ALD1_max = 0;
ALD2_max = 0;
ALD3_max = 0;
max_snow_depth = 0;

        
for (j=1:1:params_gipl.ndepth)

    if (totalPeatDepth >= params_gipl.soilLayerDepth(j))            % soil layer/node is above basal peat
        hpm_depth_index(j) = find(depth > params_gipl.soilLayerDepth(j), 1);
    else
        hpm_depth_index(j) = -99;
    end

    if (j == 1) 
        if (hpm_depth_index(j) > 0)  % peat
            mean_ann_wfps(j) = mean(annWFPS(1:hpm_depth_index(j)));
            fporos(j) = mean(poros(1:hpm_depth_index(j)));
            fpeat(j) = 1 - fporos(j);
            fmin(j) = 0;
            fwat(j) = mean(poros(1:hpm_depth_index(j))) .* annWFPS(1:hpm_depth_index(j));
            fair(j) = fporos(j) - fwat(j);
        else                         % sub-peat mineral soil
            fporos(j) = 0.65;
            fpeat(j) = 0;
            fmin(j) = 1 - fporos(j);
            fwat(j) = 0.91 * fporos(j);    %  ORIGINALLY WAS 0.71 * fporos(j);
            fair(j) = fporos(j) - fwat(j);
        end

    else
        if (hpm_depth_index(j) > 0)  % peat
            mean_ann_wfps(j) = mean(annWFPS(hpm_depth_index(j-1):hpm_depth_index(j)));
            fporos(j) = mean(poros(hpm_depth_index(j-1):hpm_depth_index(j)));
            fpeat(j) = 1 - fporos(j);
            fmin(j) = 0;
            fwat(j) = mean(poros(hpm_depth_index(j-1):hpm_depth_index(j)) .* ...
                       annWFPS(hpm_depth_index(j-1):hpm_depth_index(j)));
            fair(j) = fporos(j) - fwat(j);
        else                         % sub-peat mineral soil
            fporos(j) = 0.65;
            fpeat(j) = 0;
            fmin(j) = 1 - fporos(j);
            fwat(j) =  0.71 * fporos(j);
            fair(j) = fporos(j) - fwat(j);
        end
    end

    if (hpm_depth_index(j) > 0)  % peat
        Cond_Th(j) = 0.04 + 0.005 * mean_ann_wfps(j);   % from O?Donnel et al. 2009, via Wisser et al. 2011
        Cond_Fr(j) = 0.0141 + 0.0055 * mean_ann_wfps(j);
        Cvol_Th(j) = fpeat(j) * params_gipl.Cpeat + fwat(j) * params_gipl.Cwat;
        Cvol_Fr(j) = fpeat(j) * params_gipl.Cpeat + fwat(j) * params_gipl.Cice;  % UNFROZEN WATER CONTENT IN PEAT IS LOW (Wisser et al. 2011)
        Cvol_Sl(j) = 0.5 * (Cvol_Th(j) + Cvol_Fr(j)) + params_gipl.apparent_heat_cap * fwat(j);  % during phase change
    else                         % sub-peat mineral soil
        Cond_Th(j) = fmin(j) * params_gipl.Lmin + fwat(j) * params_gipl.Lwat + fair(j) * params_gipl.Lair;
        Cond_Fr(j) = fmin(j) * params_gipl.Lmin + fwat(j) * params_gipl.Lice + fair(j) * params_gipl.Lair;
        Cvol_Th(j) = fmin(j) * params_gipl.Cmin + fwat(j) * params_gipl.Cwat;
        Cvol_Fr(j) = fmin(j) * params_gipl.Cmin + (fwat(j) - params_gipl.Wunf) * params_gipl.Cice + params_gipl.Wunf * params_gipl.Cwat;  % UNFROZEN WATER CONTENT (= f(temp) in Wisser et al. 2011)
        Cvol_Sl(j) = 0.5 * (Cvol_Th(j) + Cvol_Fr(j)) + params_gipl.apparent_heat_cap * fwat(j);  % during phase change
    end

end

for j = params_gipl.ndepth+1:1:params_gipl.NumberOfSoilComputationNodes   % bedrock values
    if (params_gipl.soilLayerDepth(j) <= 30. )
        Cvol_Fr(j) = 1.8e6;
        Cvol_Th(j) = 1.8e6;
        fwater(j) = 0.2;
        Cvol_Sl(j) = 0.5 * (Cvol_Th(j) + Cvol_Fr(j)) + params_gipl.apparent_heat_cap * fwat(j);
        Cond_Th(j) = 2.12;
        Cond_Fr(j) = 2.54;
%        fair(j) = 0;
    else
        Cvol_Fr(j) = 2.7e6;
        Cvol_Th(j) = 2.7e6;
        fwater(j) = 0.1;
        Cvol_Sl(j) = 0.5 * (Cvol_Th(j) + Cvol_Fr(j)) + params_gipl.apparent_heat_cap * fwat(j);
        Cond_Th(j) = 2.16;
        Cond_Fr(j) = 2.51;
%        fair(j) = 0;
    end
end

% if (iyear == 3550)
%     pause;
% end

% **************************************
%  CALCULATED DAILY SNOW DEPTH, DENSITY, AND ALPHA VALUE
%     AND SAVE FOR DAILY TIME-STEP GIPL2 CALCULATIONS
% **************************************

for iday = 1:1:366
        
    T_air(iday) = daily_air_temp_time_series_linear(iday);  % °C
    P(iday) = daily_precip_time_series_linear(iday);  % mm/mo to m/mo
    snowfall(iday) = P(iday) * (T_air(iday) < -1);%  based on Willmott et al. 1995
    rainfall(iday) = P(iday) - snowfall(iday);

%        T_air(imonth) = monthly_weather_time_series.data(iyear,imonth+1);  % °C
%        P(imonth) = monthly_weather_time_series.data(iyear,imonth+13)/1000;  % mm/mo to m/mo
%        snowfall(imonth) = P(imonth) * (T_air(imonth) < -1);%  based on Willmott et al. 1995
%        rainfall(imonth) = P(imonth) - snowfall(imonth); 

%        if (T_air_month < -1) %  based on Willmott et al. 1995
%            snowfall_month = P_month;
%           rainfall_month = 0;
%        else
%            snowfall_month = 0;
%            rainfall_month = P_month;
%        end

    snowsublimation(iday) = 0.1 * swe_day_prev / 30.5;  % ???
    if (T_air(iday) > 1)
        snowmelt(iday) = (2.63 + 2.55*T_air(iday) + 0.0192*T_air(iday)*P(iday)) / 1000; % Willmott et al. 1995
    else
        snowmelt(iday) = 0;
    end
    swe(iday) = max(0, swe_day_prev + snowfall(iday) - snowmelt(iday) - snowsublimation(iday));
    swe_day_prev = swe(iday);
    imonth= min(12,floor(iday/30.5)+1);
    if (swe(iday) > 0)
        snowDepth(iday) = swe(iday) * 1000 / params_gipl.snowDensity(imonth); % 1 m water = 1000 kg
        ALFA(iday) = 1.0/params_gipl.ALFA0 + snowDepth(iday)/(0.018 + 0.00087 * params_gipl.snowDensity(imonth));
    else
        snowDepth(iday) = 0;
        ALFA(iday) = params_gipl.ALFA0;
    end
        
%        snowsublimation(imonth) = 0.1 * swe_month_prev;  % ???
%        if (T_air(imonth) > 1)
%            snowmelt(imonth) = 30.5 * (2.63 + 2.55*T_air(imonth) + 0.0192*T_air(imonth)*P(imonth)) / 1000; % Willmott et al. 1995
%        else
%            snowmelt(imonth) = 0;
%        end
%        swe(imonth) = max(0, swe_month_prev + snowfall(imonth) - snowmelt(imonth) - snowsublimation(imonth));
%        swe_month_prev = swe(imonth);
%        if (swe(imonth) > 0)
%            snowDepth(imonth) = swe(imonth) * 1000 / snowDensity(imonth); % 1 m water = 1000 kg
%            ALFA(imonth) = 1.0/ALFA0 + snowDepth(imonth)/(0.018 + 0.00087 * snowDensity(imonth));
%        else
%            snowDepth(imonth) = 0;
%            ALFA(imonth) = ALFA0;
%        end

    max_snow_depth = max(max_snow_depth,snowDepth(iday));

end  % loop through days of one year for snowpack calculations

% call/implement GIPL2 model

% ********************************************
% LOOP THROUGH 21 MONTHS OF GIPL2 SIMULATION
% ********************************************
 
for iday = 1:1:366
    
    maxABS = params_gipl.max_ABS;    % numerical iteration threshold

    soilTemp = zeros(1, params_gipl.NumberOfSoilComputationNodes);
    U1 = zeros(1, params_gipl.NumberOfSoilComputationNodes);
%        snowDepth(t)=0; edit 1-26-10
   
% 	within a loop --- time step = 1 (initial conditions)
    if  (iday == 1) % START WITH INITIAL CONDITIONS
            
        prevSoilTemp = T_init_gipl;
        U1 = prevSoilTemp;
   
    else   % not initial timestep

        U1 = prevSoilTemp;

% ---------------   
        S = 24*60*60;      %//   ! 1 day time step in seconds
%            S = 30.5*24*60*60; %//   ! 1 month (30.5 day) time step in seconds
        iter = 0; 
        
        while (iter < params_gipl.iter0 && maxABS > params_gipl.E0) %  (iter0 = 21, MAXIMUM?)
            
% 	---- computation of surface boundary coefficients  G1,G2
%            L0 = soilThermalConductivity(soilLayerDepth(1),U1(1));
%            L1 = soilThermalConductivity(soilLayerDepth(2),U1(2)); 	 
            if (U1(1) < (params_gipl.Tfr - params_gipl.FIT))
                L0 = Cond_Fr(1);
            elseif (U1(1) > (params_gipl.Tfr + params_gipl.FIT))
                L0 = Cond_Th(1);
            else
                L0 = 0.5 * (Cond_Th(1) + Cond_Fr(1));
            end
            if (U1(2) < (params_gipl.Tfr - params_gipl.FIT))
                L1 = Cond_Fr(2);
            elseif (U1(2) > (params_gipl.Tfr + params_gipl.FIT))
                L1 = Cond_Th(2);
            else
                L1 = 0.5 * (Cond_Th(2) + Cond_Fr(2));
            end

            H0 = params_gipl.soilLayerDepth(2) - params_gipl.soilLayerDepth(1);

            if snowDepth(iday) < params_gipl.E0 || prevSoilTemp(1) > params_gipl.E0 % no snow OR thawed surface%
%                if snowDepth(imonth) < E0 || prevSoilTemp(1) > E0 % no snow OR thawed surface
                G1 = 0.0;
                G2 = daily_air_temp_time_series_linear(iday);
%                G2 = monthly_weather_time_series.data(iyear,imonth+1);
                  
            elseif prevSoilTemp(1) <= 0.0 &&  snowDepth(imonth) < params_gipl.E0 % no snow and frozen surface
                G1 = 0.0;
                G2 = daily_air_temp_time_series_linear(iday);
%                G2 = monthly_weather_time_series.data(iyear,imonth+1);
 	                      
            else   % snowpack exists
    
%                ALFA = snowProperties(snowDensity, snowDepth(t));
                ALPHA = 1 / ALFA(iday);
%                ALPHA = 1 / ALFA(imonth);
%                C1 = heatCapacityDynWT(soilLayerDepth(1),prevSoilTemp(1));
                if (prevSoilTemp(1) < (params_gipl.Tfr - params_gipl.FIT))
                    C1 = Cvol_Fr(1);
                elseif (prevSoilTemp(1) > (params_gipl.Tfr + params_gipl.FIT))
                    C1 = Cvol_Th(1);
                else
                    C1 = Cvol_Sl(1);
                end
                W1 = 0.5 * (L0+L1);
                W2 = H0 * ALPHA / W1;
                W1 = 0.5 * power(H0,2) * C1 / W1 / S;
                G1 = 1.0 + W1 + W2;
                G2 = (W2 * daily_air_temp_time_series_linear(iday) + W1 * prevSoilTemp(1)) / G1;
%                G2 = (W2 * monthly_weather_time_series.data(iyear,imonth+1) + W1 * prevSoilTemp(1)) / G1;
                G1 = 1 / G1;
            end    

% 	----- Permutation and forward elimination
            P1(2) = G1;
            Q1(2) = G2;
                   
            for i = 2:1:(params_gipl.NumberOfSoilComputationNodes-1)
%                C1 = heatCapacityDynWT(soilLayerDepth(i),prevSoilTemp(i));   
                if (prevSoilTemp(i) < (params_gipl.Tfr - params_gipl.FIT))
                    C1 = Cvol_Fr(i);
                elseif (prevSoilTemp(i) > (params_gipl.Tfr + params_gipl.FIT))
                    C1 = Cvol_Th(i);
                else
                    C1 = Cvol_Sl(i);
                end
%                L2 = soilThermalConductivity(soilLayerDepth(i+1),prevSoilTemp(i+1));
                L2 = Cond_Fr(i+1) * (U1(i+1) < (params_gipl.Tfr - params_gipl.FIT)) + ... 
                      Cond_Th(i+1) * (U1(i+1) > (params_gipl.Tfr + params_gipl.FIT)) +  ...
                      0.5 * (Cond_Th(i+1) + Cond_Fr(i+1)) * ... 
                      ((U1(i+1) >= (params_gipl.Tfr - params_gipl.FIT)) && (U1(i+1) <= (params_gipl.Tfr + params_gipl.FIT)));

                H1 = params_gipl.soilLayerDepth(i+1) - params_gipl.soilLayerDepth(i);
                H2 = 0.5 * (H0 + H1);
                A1 = 0.5 * (L0 + L1) * S / C1 / (H0 * H2);
                B1 = 0.5 * (L1 + L2) * S / C1 / (H1 * H2);
                C0 = 1.0 + A1 + B1;
                P1(i+1) = B1 / (C0 - A1 * P1(i));
                Q1(i+1) = (A1 * Q1(i) + prevSoilTemp(i)) * P1(i+1) / B1;
                H0 = H1;
                L0 = L1;
                L1 = L2;
                      
%                heatCapacityOut(i,t) = C1;
%                thermalConductivityOut(i,t) = L2;
            end
%                      
% 	---- computation of the Lower boundary koef. G3 & G4
%            C1 = heatCapacityDynWT(soilLayerDepth(NumberOfSoilComputationNodes),...
%                     prevSoilTemp(NumberOfSoilComputationNodes));
            if (prevSoilTemp(params_gipl.NumberOfSoilComputationNodes) < (params_gipl.Tfr - params_gipl.FIT))
                C1 = Cvol_Fr(params_gipl.NumberOfSoilComputationNodes);
            elseif (prevSoilTemp(params_gipl.NumberOfSoilComputationNodes) > (params_gipl.Tfr + params_gipl.FIT))
                C1 = Cvol_Th(params_gipl.NumberOfSoilComputationNodes);
            else
                C1 = Cvol_Sl(params_gipl.NumberOfSoilComputationNodes);
            end
            G3 = 0.5 * power(H1,2) * C1 / L2 / S ;
            G4 = H1 * params_gipl.G0 + G3 * prevSoilTemp(params_gipl.NumberOfSoilComputationNodes);
            G3 = 1.0 / (1.0 + G3);
            G4 = G4 * G3;
                    
%  	---- Temperature computation in the last (deepest) grid node
            W1 = (G3 * Q1(params_gipl.NumberOfSoilComputationNodes) + G4) / ...
                   (1.0 - G3 * P1(params_gipl.NumberOfSoilComputationNodes));
% 	                  
            maxABS = abs(W1-U1(params_gipl.NumberOfSoilComputationNodes));
            U1(params_gipl.NumberOfSoilComputationNodes) = W1;
% 	                  
% 	 ---- Back substitution

            i = (params_gipl.NumberOfSoilComputationNodes-1);
            while (i>=1)
                W1 = P1(i+1) * U1(i+1) + Q1(i+1);
               
%	   ! check for the iterative convergence
                if (abs(W1-U1(i)) > maxABS)
                    maxABS = abs(W1-U1(i));
                end
                U1(i) = W1;
                           
                i = i - 1;

            end % END BACK-SUBSTITUTION WHILE (i>=1)                   
                       
            iter = iter + 1;
                    
        end % end while ((ITER < ITER0).AND.(maxABS > E0))
                            
        soilTemp = U1 ;

%         for (i=1:1:NumberOfSoilComputationNodes)%{//do i=1,N
        prevSoilTemp = soilTemp;                

        daily_soil_temp_out(1:params_gipl.NumberOfSoilComputationNodes,iday) = soilTemp;
%            soil_temp_out(1:NumberOfSoilComputationNodes, ((iyear - 1) * 12 + imonth)) = soilTemp';
%            soilT1 = soilTemp(1);
            
    end  % end of if loop for initial time step (if t == 1)
             
% --- ALD NOW COMPUTED ON MONTHLY MEAN TEMPS (BELOW) - MAY BE MORE STABLE (BUT WHY?)
%  ---- locate top of permafrost using 3 criteria: Tfr + FIT; Tfr; Tfr - FIT (FIT >0°C, Tfr ?0°C)
%  ----   search only to 'ndepth' soil layer, which is at 10 m; below that is considered bedrock

%     index1 = find(soilTemp(1:params_gipl.ndepth) >= (params_gipl.Tfr + params_gipl.FIT), 1, 'last');  % warm end of thawing zone (above bedrock)
%     if (isempty(index1)) 
%         ALD1_day = 0; 
%     elseif (index1 == params_gipl.ndepth)
%         ALD1_day = params_gipl.soilLayerDepth(params_gipl.ndepth);
%     else
%         ALD1_day = params_gipl.soilLayerDepth(index1) + (params_gipl.soilLayerDepth(index1+1) - params_gipl.soilLayerDepth(index1)) * ...
%              abs((params_gipl.Tfr + params_gipl.FIT) - soilTemp(index1)) / abs(eps + soilTemp(index1) - soilTemp(index1+1));
%     end
%             
%     index2 = find(soilTemp(1:params_gipl.ndepth) >= (params_gipl.Tfr), 1, 'last');  % middle of thawing zone (above bedrock)
%     if (isempty(index2)) 
%         ALD2_day = 0; 
%     elseif (index2 == params_gipl.ndepth)
%         ALD2_day = params_gipl.soilLayerDepth(params_gipl.ndepth);
%     else
%         ALD2_day = params_gipl.soilLayerDepth(index2) + (params_gipl.soilLayerDepth(index2+1) - params_gipl.soilLayerDepth(index2)) * ...
%              abs((params_gipl.Tfr) - soilTemp(index2)) / abs(eps + soilTemp(index2) - soilTemp(index2+1));
%     end
%             
%     index3 = find(soilTemp(1:params_gipl.ndepth) >= (params_gipl.Tfr - params_gipl.FIT), 1, 'last');  % cool end of thawing zone (above bedrock)
%     if (isempty(index3)) 
%         ALD3_day = 0; 
%     elseif (index3 == params_gipl.ndepth)
%         ALD3_day = params_gipl.soilLayerDepth(params_gipl.ndepth);
%     else
%         ALD3_day = params_gipl.soilLayerDepth(index3) + (params_gipl.soilLayerDepth(index3+1) - params_gipl.soilLayerDepth(index3)) * ...
%              abs((params_gipl.Tfr - params_gipl.FIT) - soilTemp(index3)) / abs(eps + soilTemp(index3) - soilTemp(index3+1));
%     end
% 
%     ALD1_max = max(ALD1_day, ALD1_max);
%     ALD2_max = max(ALD2_day, ALD2_max);
%     ALD3_max = max(ALD3_day, ALD3_max);

%  ---- Save daily temps for computing monthly values at end of year

    daily_soil_temps(:,iday) = soilTemp';
    
end  % loop of for loop through days (for iday = 1:366)

% ALD1 = ALD1_max;
% ALD2 = ALD2_max;
% ALD3 = ALD3_max;

%  ---- soil layer temps equal average of the upper and lower nodal temps

daily_layer_temps = 0.5*(daily_soil_temps(1:params_gipl.ndepth-1,:) + daily_soil_temps(2:params_gipl.ndepth,:));
  
for jmonth=1:12   % 30 day months, starting with DOY 4
    monthly_layer_temps(:,jmonth) = mean(daily_layer_temps(:,(jmonth-1)*30+4:jmonth*30+3),2);
    monthly_soil_temps(:,jmonth) = mean(daily_soil_temps(:,(jmonth-1)*30+4:jmonth*30+3),2);
end

soil_layer_temp_month = monthly_layer_temps';
soil_node_temp_month = monthly_soil_temps';

% --- AT END OF YEAR, COMPUTE MAX ALDs FROM MONTHLY SOIL NODE TEMPS

ALD1_max = 0;
ALD2_max = 0;
ALD3_max = 0;
ALD1_month = 0;
ALD2_month = 0;
ALD3_month = 0;
monthly_node_T = zeros(params_gipl.ndepth,1);

for jmonth=1:12
    monthly_node_T = monthly_soil_temps(1:params_gipl.ndepth,jmonth);
    if (monthly_soil_temps(1,jmonth) > params_gipl.Tfr + params_gipl.FIT)
        index1 = find(monthly_node_T > (params_gipl.Tfr + params_gipl.FIT), 1, 'last');  % warm end of thawing zone (above bedrock)
        if (isempty(index1)) 
            ALD1_month = 0; 
        elseif (index1 == params_gipl.ndepth)
            ALD1_month = params_gipl.soilLayerDepth(params_gipl.ndepth);
        else
            ALD1_month = params_gipl.soilLayerDepth(index1) + (params_gipl.soilLayerDepth(index1+1) - params_gipl.soilLayerDepth(index1)) * ...
                 abs((params_gipl.Tfr + params_gipl.FIT) - monthly_node_T(index1)) / abs(eps + monthly_node_T(index1) - monthly_node_T(index1+1));
        end
    else
        ALD1_month = 0;
    end

    if (monthly_soil_temps(1,jmonth) > params_gipl.Tfr)
        index2 = find(monthly_node_T > (params_gipl.Tfr), 1, 'last');  % middle of thawing zone (above bedrock)
        if (isempty(index2)) 
            ALD2_month = 0; 
        elseif (index2 == params_gipl.ndepth)
            ALD2_month = params_gipl.soilLayerDepth(params_gipl.ndepth);
        else
            ALD2_month = params_gipl.soilLayerDepth(index2) + (params_gipl.soilLayerDepth(index2+1) - params_gipl.soilLayerDepth(index2)) * ...
                 abs((params_gipl.Tfr) - monthly_node_T(index2)) / abs(eps + monthly_node_T(index2) - monthly_node_T(index2+1));
        end
    else
        ALD2_month = 0;
    end

    if (monthly_soil_temps(1,jmonth) > params_gipl.Tfr - params_gipl.FIT)
        index3 = find(monthly_node_T > (params_gipl.Tfr - params_gipl.FIT), 1, 'last');  % cool end of thawing zone (above bedrock)
        if (isempty(index3)) 
            ALD3_month = 0; 
        elseif (index3 == params_gipl.ndepth)
            ALD3_month = params_gipl.soilLayerDepth(params_gipl.ndepth);
        else
            ALD3_month = params_gipl.soilLayerDepth(index3) + (params_gipl.soilLayerDepth(index3+1) - params_gipl.soilLayerDepth(index3)) * ...
                 abs((params_gipl.Tfr - params_gipl.FIT) - monthly_node_T(index3)) / abs(eps + monthly_node_T(index3) - monthly_node_T(index3+1));
        end
    else
        ALD3_month = 0;
    end
    
    ALD1_max = max(ALD1_month, ALD1_max);
    ALD2_max = max(ALD2_month, ALD2_max);
    ALD3_max = max(ALD3_month, ALD3_max);
end

ALD1 = ALD1_max;
ALD2 = ALD2_max;
ALD3 = ALD3_max;

% --- assemble soil temp output into array year, month, soil_level_above_bedrock

% monthly_soil_temp_array = zeros(num_years,12,ndepth-1);
% monthly_soil_temp_series = zeros(ndepth-1,num_years*12);
% for jyear = 1:1:num_years
%     for jmonth = 1:1:12
%         for jlevel = 1:1:ndepth-1
%             monthly_soil_temp_array(jyear,jmonth,jlevel) = mean(mean(soil_temp_out(jlevel:jlevel+1,((jyear-1)*366+(jmonth-1)*30+1):((jyear-1)*366+(jmonth)*30))));
%             monthly_soil_temp_series(jlevel,(jyear-1)*12+jmonth) = mean(mean(soil_temp_out(jlevel:jlevel+1,((jyear-1)*366+(jmonth-1)*30+1):((jyear-1)*366+(jmonth)*30))));
%         end
%     end
% end

        