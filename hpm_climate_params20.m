%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% READ IN MONTHLY WEATHER DATA
%    ? for other sites, monthly file generated externally
%   also can generate monthly climate data using sinusoidal curves plus
%   white noise.
%   main code requires monthly climate data truncated to the simulation
%   period. This processing should be done in this file.
%   Simulation period is in y BP, where 0 is 1950 CE and -50 is 2000 CE.
%   outputs are precip_forcing_mon, temp_forcing_mon
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%-------------------------------
%    Set annual temperature & precip to use in ET calcs. 
%-------------------------------
ann_temp = -10; % site mean annual temp (C)
ann_ppt = 0.35;    % long-term mean precipitation (m/y)

% Need to ask Steve what these are about? The only place they seem to be
% used is in the param output file. Setting to 0.
ann_temp_orig = ann_temp;
ann_temp_offset = 0; %+2;
ann_temp = ann_temp + ann_temp_offset;

ann_ppt_orig = ann_ppt;
% ann_ppt = ann_ppt * exp(0.068 *(ann_temp - ann_temp_orig));
mon_ppt_offset = 0; %0.060; % m per MONTH!
annppt_modifier = 0.;   % used in main code, set a default value.

%-------------------------------
%    Read in monthly climate data 
%-------------------------------
climForcingFile = 'toolik_monthly_T_P.csv'; 
% climForcingFile = name with specifics if making climate file
T = readtable(climForcingFile);
TA = table2array(T);

years_forcing = TA(:,1);
months_forcing = TA(:,2);

% Get years_forcing into common format: oldest to youngest years, in BP,
% where 0 BP = 1950 CE and -50 BP = 2000 CE
% years_BP_forcing = -years_forcing; % 0 BP = 1950 CE 
years_BP_forcing = -years_forcing + 4001; % 0 BP = 1950 CE 

%-------------------------------
%    TRUNCATE CLIMATE DATA TO SIMULATION PERIOD
%-------------------------------

index_forcing_start = find(years_BP_forcing <= sim_start,1);
index_forcing_end = find(years_BP_forcing <= sim_end,1) + 11; % 12 months per year in file

years_forcing_ann = years_forcing(index_forcing_start:index_forcing_end);
months_forcing_ann = months_forcing(index_forcing_start:index_forcing_end);

precip_forcing_mon = TA(index_forcing_start:index_forcing_end,3); % in m;
temp_forcing_mon = TA(index_forcing_start:index_forcing_end,5); % in C;

% FOR MER BLEUE GFDL DRIVERS FROM ZACK
% if (params.gfdl_model_flag < 1.5)
%     precip_forcing_mon = precipforc_CM3(index_forcing_start:index_forcing_end) / 1000.;  % mm/mon to m/mon
%     temp_forcing_mon = tsaforc_CM3(index_forcing_start:index_forcing_end) - 273.15;      % deg K to °C
%  else
%     precip_forcing_mon = precipforc_ESM2M(index_forcing_start:index_forcing_end) / 1000.;  % mm/mon to m/mon
%     temp_forcing_mon = tsaforc_ESM2M(index_forcing_start:index_forcing_end) - 273.15;      % deg K to °C
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permafrost thaw/refreeze ramps
% del_T_thaw_freeze = zeros(49812,1);
% index_vector = 1:1:49812;
% del_T_thaw_freeze(9001:24000) = (index_vector(9001:24000) - index_vector(9001)) / ...
%     (index_vector(24000) - index_vector(9001)) * 15.;
% del_T_thaw_freeze(24001:39000) = (1 - abs((index_vector(24000) - index_vector(24001:39000)) / ...
%     (index_vector(24000) - index_vector(39000)))) * 15.;
% 
% temp_forcing_mon = TA(:,5) + del_T_thaw_freeze;
% end permafrost thaw/refreeze ramps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permafrost thaw/refreeze small amp sinusoidal around ~0°C
% del_T_thaw_freeze = zeros(49812,1);
% index_vector = 1:1:49812;
% del_T_thaw_freeze(9001:15000) = (index_vector(9001:15000) - index_vector(9001)) / ...
%      (index_vector(15000) - index_vector(9001)) * 10.;
% del_T_thaw_freeze(15001:33000) = 10. + sin((index_vector(15001:33000)-9001)/6000 * 2*pi) * 3.;
% del_T_thaw_freeze(33001:39000) = (1 - abs((index_vector(33000) - index_vector(33001:39000)) / ...
%      (index_vector(33000) - index_vector(39000)))) * 10.;
% 
% temp_forcing_mon = TA(:,5) + del_T_thaw_freeze;
% end permafrost thaw/refreeze ramps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% permafrost thaw/refreeze large sinusoidal
% del_T_thaw_freeze = zeros(49812,1);
% index_vector = 1:1:49812;
% del_T_thaw_freeze(9001:39000) = sin((index_vector(9001:39000)-9001)/20000 * 2*pi).^2 * 15.;
% 
% temp_forcing_mon = TA(:,5) + del_T_thaw_freeze;
% end permafrost thaw/refreeze ramps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------
%    Use to generate synthetic temperature and precip data.
%    ann_temp_amp is saved to params.x; add other values back into
%    parameter file from bottom if using synthetic precip, as well as
%    adding fprintf back to main code.
%-------------------------------
ann_temp_amp = 0; % site monthly temp amplitude (half max minus min)
%    % this is used to generate a sinusoidal monthly air temperature
%     for jmonth = 1:12
%         temp_forcing_mon(jmonth) = T_mean + ann_temp_amp * sin(2*pi*jmonth/12 - 2*pi/3);
%     end
%    precip_forcing_mon = ann_ppt/12;  % to be replaced with monthly input data from RENATO

% % these 'temp_xx' parameters are used to generate noise in holocene temp
% temp_flag = 14;   % permafrost case to generate synthetic temps
% temp_rand_persist = 0.5; %0.95;   %random persistence parameter
% temp_amp1 = 2.;   % max. amplitude of random noise in mean annual temp
% temp_amp2 = 3.;   % max. amplitude of random noise in annual temp range
%    
% 
% ppt_flag = 14;   % 1 - sinusoidal (or const if amps are 0); 3 - ramp; 5 - ramp up/down or down/up; 9 - MB reconstruction; 11 = abrupt drying, then gradual rewetting
%                 % 13 - Ner Bleue serge muller reconstruction; 14 = permafrost
%                 
% ppt_amp1 = 0.133;  % amplitude (m/y) of sine or ramp (+: ramp up; -: ramp down) (+: ramp up then down; -: ramp down then up)
% ppt_amp2 = 0.0;  % random noise amplitude (m/y) (values ~0.1-ish seem resaonable) 
% % NOTE: for Muller MB reconstruction amp1 ~ 2 is variability multiplier to approx match Muller confidence intervale)
% % note: for MB ppt reconstruction, ppt_amp1 is noise during constant ppt
% %       intervals and ppt_amp2 is noise during ppt transitions (zicheng thinks amp2 > amp1)
% 
% ppt_rand_persist = 0.94;   % persistence of random ppt variability (see main code just before loop through yearss)
% % NOTE: for Muller MB reconstruction rand_persist ~ 0.97 is variability factor to approx match Muller confidence intervale)

%%%%%%%%%%%%% Add these to the run parameters output section in the main code, add back into params file if using synthetic temperature. 
% fprintf(fid4,'1-sine,3-ramp,5-ramps,9-MB %6.2f \n',params.ppt_flag);
% fprintf(fid4,'sine/ramp amp [m/y]        %6.2f \n',params.ppt_amp1);
% fprintf(fid4,'ppt noise amp [m/y]        %6.2f \n',params.ppt_amp2);
% fprintf(fid4,'ppt_noise_persist          %6.3f \n',params.ppt_rand_persist);


% --- XXX bump up temp by 6 degrees for a quick test  
%    air_temp_month = air_temp_month + 6;
    
%    gdd_sum(itime) = sum(air_temp_month > 0) * 30.5; % degree-days (NOT USED)
    
% Moved from main code 
% [REC_anntemp REC_annppt] = hpm_climate10(istep, params, num_years);  %reconstructed annual climate

%REC_annppt = importdata('HPM_CAR_LPX_precip.csv');
% REC_annppt = importdata('HPM_CAR_t1g_precip.csv');
%  REC_annppt = REC_annppt / 1000. ; % convert mm/y to m/y  NOT NECESSARY


% Mer Bleue reconstruction BUT 10°C COLDER and 40% LESS PRECIPITATION!!!!
%    REC_anntemp = REC_anntemp + params.ann_temp_offset;
%    REC_annppt = REC_annppt * exp(0.068 *(params.ann_temp - params.ann_temp_orig));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------
% accumulate monthly precip to annual
%-------------------------------
jan = 1:12:12*sim_len;
feb = 2:12:12*sim_len;
mar = 3:12:12*sim_len;
apr = 4:12:12*sim_len;
may = 5:12:12*sim_len;
jun = 6:12:12*sim_len;
jul = 7:12:12*sim_len;
aug = 8:12:12*sim_len;
sep = 9:12:12*sim_len;
oct = 10:12:12*sim_len;
nov = 11:12:12*sim_len;
dec = 12:12:12*sim_len;

precip_forcing_ann = precip_forcing_mon(jan) + precip_forcing_mon(feb) + ...
    precip_forcing_mon(mar) + precip_forcing_mon(apr) + precip_forcing_mon(may) + ...
    precip_forcing_mon(jun) + precip_forcing_mon(jul) + precip_forcing_mon(aug) + ...
    precip_forcing_mon(sep) + precip_forcing_mon(oct) + ...
    precip_forcing_mon(nov) + precip_forcing_mon(dec);

temp_forcing_ann = (temp_forcing_mon(jan) + temp_forcing_mon(feb) + ...
    temp_forcing_mon(mar) + temp_forcing_mon(apr) + temp_forcing_mon(may) + ...
    temp_forcing_mon(jun) + temp_forcing_mon(jul) + temp_forcing_mon(aug) + ...
    temp_forcing_mon(sep) + temp_forcing_mon(oct) + ...
    temp_forcing_mon(nov) + temp_forcing_mon(dec)) / 12.;
