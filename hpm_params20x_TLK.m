% Parameter values for HPM20;  S. Frolking, J. Talbot, 2008-2015

% v.20: making version for version control - June 2015

% v.13: re-organized to make changing # of PFTs easier (Sept. 2014); set code to take arbitrary number
% of PFTs.  Added old/new Carbon (double PFTs).  

% v. 12: modified to read in atmosphere del-14C and more flexible simulation times; summer 2014

% v.8: adjusted lawn and hummock sphagnum NPP sensitivity to WTD

% new in version 6: WTD_range and PD_range have 2 sets of values (see 'vegetation_NPP_4.m')

% **************
%  Output file base name
% **************

outname = 'HPM-20 output files\HPM_20_soiltemp_TLK_small_SINE_CCSM4_RCP85_to_2100_NPPQ10_15_5_ETtest';

sim_start = 500; % years BP (before 'present'), where 0 PB = 1950 CE
sim_end = -1;   % years BP  (-150 BP = 2100 CE)
% sim_start = 8500; % years BP (before 'present'), where 0 PB = 1950 CE
% sim_end = -550;   % years BP  (-550 BP = 2500 CE)
sim_len = sim_start - sim_end + 1;  % simulation length (years)

gipl_flag = 1; % if 0 (or 1) skip (or run) GIPL soil physics model: no (or yes) temperature effect on decomp
pf_flag = 1; % if 1 site has permafrost; otherwise 0 

% **************
%  SITE CLIMATE
% **************
% read in monthly climate data from climate processing file
hpm_climate_params20;

% gfdl_model_flag = 2;  % FOR MER BLEUE RUNS: if 1: GFDL CM3 RCP8.5; if 2: GFDL ESM2M RCP8.5
gfdl_model_flag = 1;  % FOR TOOLIK RUNS: if 1: CCSM4 RCP8.5; if 2: CCSM4 RCP4.5

flag_RCP = 1; % 1 = RCP8.5, 2 = RCP4.5

ET_0 = 0.2;      % base evapotranspiration (m/y)
ET_0 = 0.281*exp(0.068 * ann_temp);

ald_0 = 0.3;  % first year active layer depth, if needed (m)
wtd_0 = 0.07; % initialization period water table depth (m)
start_depth = 0.25; % depth of initial peat accumulation (m) at which water balance calculations begin
% start_depth = 0.5; % depth of initial peat accumulation (m) at which water balance calculations begin


% **************
%  BOG or FEN??
% **************

bog_fen_id = 2; % fen-to-bog = 1, persistent fen or permafrost = 2

if (bog_fen_id <1.5)   % FEN-TO-BOG VALUES
    
%    Roff_c1 = max(0,ann_ppt - ET_0 + 0.1); % max runoff + max ET = mean annual precip + 0.1 m/y
    Roff_c1 = max(0,ann_ppt - ET_0 + 0.05); % max runoff + max ET = mean annual precip + 0.05 m/y
    Roff_c2 = 0.2;  % linear increase in runoff (m/y) per meter of total peat height
    Roff_c2a = 1.2 * start_depth;  % peat height needed to get base run-off (factor = 1 +c2*(H-c2a))

    anoxia_scale_length = 0.3;  % exponential decline in decomp in catotelm from wfps_sat_rate to wfps_min_rate
%    anoxia_scale_length = 0.25;  % exponential decline in decomp in catotelm from wfps_sat_rate to wfps_min_rate

    runon_c1 = 1.0;  % total peat depth (m) where runon declines by ~50%
    runon_c2 = 0.5;  % controls rate of decline of runon (see 'HPM vegetation productivity.xls')
    runon_c3 = 0.0;   % magnitude of maximum runon (m/y)

% else   %  PERENNIAL FEN VALUES
%     
% 
%     Roff_c1 = max(0,ann_ppt - ET_0 + 0.1); % max runoff + max ET = mean annual precip
%     Roff_c2 = 0.02;  % linear increase in runoff (m/y) per meter of total peat height
%     Roff_c2a = 1.;  % peat height needed to get base run-off (factor = 1 +c2*(H-c2a))
% 
%     anoxia_scale_length = 3.0;  % exponential decline in decomp in catotelm from wfps_sat_rate to wfps_min_rate
% 
%     runon_c1 = 5.0;  % total peat depth (m) where runon declines by ~50%
%     runon_c2 = 0.5;  % controls rate of decline of runon (see 'HPM vegetation productivity.xls')
%     runon_c3 = 0.3; % magnitude of maximum runon (m/y)

else   %  PERMAFROST SITE VALUES
    

    Roff_c1 = max(0,ann_ppt - ET_0 + 0.1); % max runoff + max ET = mean annual precip
    Roff_c2 = 0.2;  % linear increase in runoff (m/y) per meter of total peat height
    Roff_c2a = 1.;  % peat height needed to get base run-off (factor = 1 +c2*(H-c2a))
    Roff_c2a = 1.2 * start_depth;  % peat height needed to get base run-off (factor = 1 +c2*(H-c2a))

    anoxia_scale_length = 0.3;  % exponential decline in decomp in catotelm from wfps_sat_rate to wfps_min_rate

    runon_c1 = 1.0;  % total peat depth (m) where runon declines by ~50%
    runon_c2 = 0.5;  % controls rate of decline of runon (see 'HPM vegetation productivity.xls')
    runon_c3 = 0.0; % magnitude of maximum runon (m/y)

end

% ???  ***SHOULD Roff_c1 EVER BE ZERO???***

% **************
%  VEGETATION 
% **************
%   NOTE: plants don't grow, so litterfall = NPP

% RUN WITH DOUBLE PFTS FOR OLD-NEW CARBON ANALYSIS

tf_old_new = 0; % 1: double PFTs for old/new; otherwise = 0 & do not do this
tf_old_new_timing = 150;  % years before end of simulation to switch 

if (tf_old_new > 0.5)
    year_old2new = sim_len - tf_old_new_timing;
else
    year_old2new = sim_len + tf_old_new_timing;  % i.e., never happens
end

% ARCTIC VERSION  
%    arctic version will use active layer depth rather than peat height for NPP

num_veg = 5;
 
PFT_1_name = 'moss'; 
PFT_2_name = 'sedge_ag'; 
PFT_3_name = 'sedge_bg'; 
PFT_4_name = 'shrub_ag'; 
PFT_5_name = 'shrub_bg'; 
 
mosses =    [ 1 0 0 0 0 ];
vasculars = [ 0 1 1 1 1 ];
sedges =    [ 0 1 1 0 0 ];
woody =     [ 0 0 0 1 1 ];

PFT_param = zeros(num_veg,12);
 
% *** PFT Parameters                   ** PD not used **
%                 WTD_0, WTD_-, WTD_+, PD_0, PD_-, PD_+, ALD_0, ALD_-, ALD_+, NPP_rel, NPP_AG, k_exp   
PFT_param(1,:) = [ 0.1   0.09    0.35   1.0   2.   19.   1.0    19.    29.     0.5      1.0     0.04  ]; % moss
PFT_param(2,:) = [ 0.025 0.15    0.20   1.0   2.   19.   1.5    1.0    29.     1.0      1.0     0.25  ]; % sedge aboveground
PFT_param(3,:) = [ 0.025 0.15    0.20   1.0   2.   19.   1.5    1.0    29.     1.0      0.0     0.225 ]; % sedge belowground
PFT_param(4,:) = [ 0.3   0.15    3.5    1.0   2.   19.   2.0    1.5    29.     1.3      1.0     0.15  ]; % shrub aboveground
PFT_param(5,:) = [ 0.3   0.15    3.5    1.0   2.   19.   2.0    1.5    29.     0.7      0.0     0.10  ]; % shrub belowground

% 'tf_xxx' parameters below are used if vascular PFTs are partitioned into two
% component PFTs (aboveground and belowground), with some different
% parameters (e.g., k_exp)

% tf_sedge_root =     [ 0 0 0 1 0 ];
% tf_non_sedge_root = [ 0 0 0 0 1 ];
% tf_moss =           [ 1 0 0 0 0 ];
% tf_ag_pft =         [ 1 1 1 1 1 ] - tf_sedge_root - tf_non_sedge_root;

% NON-ARCTIC VERSION
%    non-arctic version will use peat height rather than active layer depth for NPP

% num_veg = 13;
% 
% PFT_1_name = 'min_grass';
% PFT_2_name = 'min_herb'; 
% PFT_3_name = 'min_sedge'; 
% PFT_4_name = 'decid_shrub'; 
% PFT_5_name = 'brown_moss'; 
% PFT_6_name = 'hollow_sphagnum'; 
% PFT_7_name = 'lawn_sphagnum';
% PFT_8_name = 'hummock_sphagnum'; 
% PFT_9_name = 'feather_moss'; 
% PFT_10_name = 'omb_herb'; 
% PFT_11_name = 'omb_sedge'; 
% PFT_12_name = 'evrgn_shrub'; 
% PFT_13_name = 'tree'; 
% 
% mosses =    [ 0 0 0 0 1 1 1 1 1 0 0 0 0 ];
% vasculars = [ 1 1 1 1 0 0 0 0 0 1 1 1 1 ];
% sedges =    [ 0 0 1 0 0 0 0 0 0 0 1 0 0 ];
% woody =     [ 0 0 0 1 0 0 0 0 0 0 0 1 1 ];

%reducing PFT number to try old/new C run at MB

% num_veg = 10;
% 
% PFT_1_name = 'min_grass';
% PFT_2_name = 'min_sedge'; 
% PFT_3_name = 'decid_shrub'; 
% PFT_4_name = 'lawn_sphagnum';
% PFT_5_name = 'hummock_sphagnum'; 
% PFT_6_name = 'feather_moss'; 
% PFT_7_name = 'omb_herb'; 
% PFT_8_name = 'omb_sedge'; 
% PFT_9_name = 'evrgn_shrub'; 
% PFT_10_name = 'tree'; 
% 
% mosses =    [ 0 0 0 1 1 1 0 0 0 0 ];
% vasculars = [ 1 1 1 0 0 0 1 1 1 1 ];
% sedges =    [ 0 1 0 0 0 0 0 1 0 0 ];
% woody =     [ 0 0 1 0 0 0 0 0 1 1 ];

% PFT_param = zeros(num_veg,12);
% 
% % *** PFT Parameters                                    ***ALD NOT USED****
% %                  WTD_0, WTD_-, WTD_+, PD_0, PD_-, PD_+, ALD_0, ALD_-, ALD_+, NPP_rel, NPP_AG, k_exp   
% PFT_param(1,:) =  [ 0.4    0.4    0.4    0.01  1.0   1.0   0.1    1.     1.     3*0.75   0.5     0.2  ]; % minerotrophic grass
% PFT_param(2,:) =  [ 0.1    0.3    0.3    0.3   1.0   1.0   0.1    1.     1.     3*1.0    0.2     0.4  ]; % minerotrophic herb
% PFT_param(3,:) =  [ 0.1    0.4    0.4    0.1   2.0   2.0   0.1    1.     1.     3*1.0    0.2     0.3  ]; % minerotrophic sedge
% PFT_param(4,:) =  [ 0.2    0.2    1.0    1.0   2.0   2.0   0.1    1.     1.     3*0.5    0.5     0.25 ]; % deciduous shrub
% PFT_param(5,:) =  [ 0.01   0.2    0.05   0.1   1.5   1.5   0.1    1.     1.     1*0.5    1.0     0.1  ]; % brown moss
% PFT_param(6,:) =  [ 0.01   0.2    0.05   2.0   1.0   19.   0.1    1.     1.     1*0.5    1.0     0.1  ]; % hollow sphagnum
% PFT_param(7,:) =  [ 0.1    0.3    0.4    2.0   1.0   19.   0.1    1.     1.     1*0.5    1.0     0.07 ]; % lawn sphagnum
% PFT_param(8,:) =  [ 0.2    0.1    0.5    2.0   1.0   19.   0.1    1.     1.     1*0.5    1.0     0.05 ]; % hummock sphagnum
% PFT_param(9,:) =  [ 0.4    0.4    0.6    4.0   6.0   19.   0.1    1.     1.     1*0.25   1.0     0.1  ]; % feathermoss
% PFT_param(10,:) = [ 0.2    0.2    0.2    4.0   2.0   19.   0.1    1.     1.     1*0.25   0.5     0.3  ]; % ombrotrophic herb
% PFT_param(11,:) = [ 0.2    0.3    0.3    4.0   2.0   19.   0.1    1.     1.     1*0.5    0.2     0.2  ]; % ombrotrophic sedge
% PFT_param(12,:) = [ 0.3    0.3    1.0    4.0   2.0   19.   0.1    1.     1.     1*0.5    0.5     0.2  ]; % evergreen shrub
% PFT_param(13,:) = [ 0.8    0.3    10.0   2.0   20.   10.   0.1    1.     1.     1*2.8    0.7     0.3  ]; % tree

% **************
% ALL VERSIONS
% **************

for ii = 1:1:num_veg
    WTD_opt(ii)     = PFT_param(ii,1);      % optimal water table depth (m)
    WTD_range(1,ii) = PFT_param(ii,2);      % variance on shallow WTD side (m)
    WTD_range(2,ii) = PFT_param(ii,3);      % variance on deep WTD side (m)
    PD_opt(ii)      = PFT_param(ii,4);      % optimal peat height (m)
    PD_range(1,ii)  = PFT_param(ii,5);      % variance on shallow PD side (m)
    PD_range(2,ii)  = PFT_param(ii,6);      % variance on deep PD side (m)
    ALD_opt(ii)     = PFT_param(ii,7);      % optimal active layer depth (m)
    ALD_range(1,ii) = PFT_param(ii,8);      % variance on shallow ALD side (m)
    ALD_range(2,ii) = PFT_param(ii,9);      % variance on deep ALD side (m)
    NPP_rel(ii)     = PFT_param(ii,10);     % relative NPP
    ag_frac_npp(ii) = PFT_param(ii,11);     % aboveground fraction of NPP
    bg_frac_npp(ii) = 1. - ag_frac_npp(ii); % belowground fraction of NPP
    k_exp(ii)       = PFT_param(ii,12);     % litterbag k-value (1/year)
end

% ************************
%  DECOMP PARAMETERS
% ************************

% **************
% ARCTIC VERSION  
% **************
k_exp_temp = 4.;   % litter incubation temperature (ฐC)from Hobbie paper

k_0 = k_exp .* (1 + 3 * k_exp);  %see spreadsheet 'simple decomp models.xls'; adjusts k_0 for m/m0 model of decay
k_month_0 = k_0 / 0.411 / 12;  % convert 'per-0.4-yr' (Hobbie study) to 'per-year' to 'per-month'

% ???   does this make sense (next line)?????
% k_month_0 = k_0 / 12 * (2 + 3*exp(-ann_temp/9)).^(-(ann_temp - 10)/10);  % convert year at 6ฐ (mer bleue) to year at 10ฐ (ref)
k_month_0 = k_month_0 * (2 + 3*exp(-k_exp_temp/9)).^(-(k_exp_temp - 10)/10);  % convert month at 4ฐ (Hobbie incubation) to year at 10ฐ (ref)

% tt = 5; % years of decomp to match algorithms (5 is reasonable, i.e., long for litter bags)
% k_0 = k_exp * (1 + ((exp(k_exp*tt) - 1) ./ (k_exp * tt) -1);

% NPP_rel = NPP_rel .* [3 3 3 3 3 1 1 1 1 1 1 1 1];

% ONLY MOSSES
%NPP_rel   = NPP_rel .* mosses;

% ONLY VASCULAR
% NPP_rel   = NPP_rel .* vascular;

max_npp = 0.75 * 2.25;   % approximate absolute maximum total NPP for all vegetation at mean annual T = 10ฐC, kg/m2/y
                         % this if for TOOLIK (ann_temp = -10C=ฐC) to get final value of 0.75 kg/m2/y
% max_npp = 3.23;        % approximate absolute maximum total NPP for all vegetation at mean annual T = 10ฐC, kg/m2/y
% max_npp = 2.75;        % approximate maximum total NPP for all vegetation, kg/m2/y

q10_npp = 1.5;   % see Julie Talbot email of 4 June 2014
max_npp = max_npp * q10_npp^((ann_temp - 10)/10)

% total_npp = hpm_vegNPP10(WTD_opt,WTD_range,PD_opt,PD_range,NPP_rel,k_0);

%--- BUILD TOTAL NPP SURFACE
%      MER BLEUE USES WTD AND PD
%      TOOLIK USES WTD AND ALD

% total_npp = hpm_vegNPP12(WTD_opt,WTD_range,PD_opt,PD_range,NPP_rel,k_0,num_veg); 
total_npp = hpm_vegNPP20(WTD_opt,WTD_range,ALD_opt,ALD_range,NPP_rel,k_0,num_veg);

NPP_rel = NPP_rel * (max_npp / total_npp)   % scale relative NPP so that max sum NPP ~ 'max_npp'

NPP_Q10 = q10_npp;  % not used; using Q10(T) formulation (as about 7 lines up)

lag_years = 10;  % # years averaging WTD for vascular plants (1 year for non-vascular)

disp(sprintf('numveg: %d  total NPP (kg/m2/y): %d', num_veg, total_npp));
disp(sprintf('ann_ppt (m/y): %d  ann_ET0 (m/y): %d, ann_runoff0 (m/y): %d', ann_ppt, ET_0, Roff_c1));

% ****************************
%   DECOMPOSITION
% ****************************

% initial decomposition (mass-loss) rates and anoxia factor 
%   (make anoxia factor more variable, as in new paper by Blodau?)

wfps_opt = 0.45;  % must be <= 0.5; optimum WFPS for decomposition (see speadsheet 'simple decomp models.xls')
wfps_max_rate = 1.0;   % decomp rate multiplier at WFPS = WFPS_opt.
wfps_sat = 1.0;      % WFPS at saturation
wfps_sat_rate = 0.3;   % decomp rate multiplier at WFPS = 1.0 (i.e., at annual mean WTD).
wfps_min_rate = 0.001;   % decomp rate multiplier minimum,deep in catotelm.
wfps_curve = (wfps_sat - wfps_opt)^2 / (4 * (wfps_max_rate - wfps_sat_rate)); % parabola with value of 0.1 at WFPS = 1.0

% ****************************
% ROOT PROFILES AND LITTER INPUT
% ****************************

% FOLLOWING BAUER (2004)
%   sedge root litter input: exponential decay profile to 2 m; 80% above 0.2 to 0.3 m.
%          litter_input [kg/m2/m] = beta * exp(-alpha * depth) * layer_thickness
%   other vascular root litter: ~constant profile to max(WTD, 0.2 m)
%          litter_input [kg/m2/m] = [bg_npp / max(WTD, 0.2)] * layer_thickness
%   moss root litter input equals zero

% c3 = approx, minimum depth [m] for non-sedge root profile ((most roots above max of WTD & c3)
% c5 = parameter [--] controlling steepness of transition from roots to no roots for non-sedge profile
%    c5 used for smooth transition in non-sedge root zone, otherwise simpler abrupt transition is used.
% c4 = maximum depth [m] for sedge root profile (root profile is exponential decay to this depth)
% d80 = depth in meters to 80% of sedge roots
% alpha [1/m]; controls exponential decay in sedge root profile; denominator is 1 minus fraction of roots above 'd80'

rootin_c3 = 0.2;  
rootin_d80 = 0.3;          
rootin_alpha = -log(rootin_d80) / (1. - 0.8);   
rootin_c4 = 2.0;  
rootin_c5 = 0.04;

% **************
%   BULK DENSITY
% **************

min_bulk_dens = 50.;   % kg/m3
del_bulk_dens = 80.;   % bulk density increase down profile, kg/m3
dens_c1 = 0.333;  % m_star value at which bulk density rises halfway from min to max
dens_c2 = 0.20;  % parameter controlling steepness of bulk density transition (smaller is steeper)
OM_dens = 1300; % density of organic matter [kg/m3]

% ****************************
%  WATER BALANCE
% ****************************

% ET_wtd_1 = 0.1;   % WTD threshold for full ET (m)
% ET_wtd_2 = 0.5;   % WTD threshold for low ET (m)
% ET_param = 1/(ET_wtd_2 - ET_wtd_1);   % linear drop in ET as WTD drops from ET_wtd_1 to ET_wtd_2

% modifications June 2008 based on Lafleur et al. 2005 ET from Mer Bleue paper
ET_wtd_1 = 0.1;   % WTD threshold for full ET (m)
ET_wtd_2 = 0.5;   % WTD threshold for low ET (m)
ET_wtd_3 = 1.5;   % min ET = ET_0 รท ET_wtd_3
ET_param = 1/(ET_wtd_2 - ET_wtd_1);   % linear drop in ET as WTD drops from ET_wtd_1 to ET_wtd_2

ET_temp_sens = 0.05; % 5% change in ET per degree C (Brummer et al. 2011; Ag. For. Met. 153, 14-30)
            % replaced with ET_0(T) function from J Talbot (see spreadsheet)

% for value of Roff_c1 see top of this script.
Roff_c3 = 0.5; % minimum profile relative transmissivity (see hpm_WatBal7.m')
Roff_c4 = -0.0; % threshold water table depth for extra spillover (= Roff_c4 - WTD)

%  runon_c1 = 10.0;  % total peat depth (m) where runon declines by ~50%
% runon_c2 = 0.5;  % controls rate of decline of runon (see 'HPM vegetation productivity.xls')
% runon_c3 = 0.2; % magnitude of maximum runon (m/y)

% calculates peat cohort fractional water content (0 - 1) above the water table
%   depends on distance above water table and peat bulk density;

% Parameters (see excel spreadsheet: 'anoxia & bulk dens & WFPS % profile.xls')

wfps_c1 = 0.03;
wfps_c2 = 0.5;
wfps_c3 = 20;
% wfps_c2 = 0.7;
% wfps_c3 = 60;

if (tf_old_new > 0.5) 
    num_veg = 2. * num_veg;
    PD_opt = [PD_opt PD_opt];
    PD_range = [PD_range PD_range];
    ALD_opt = [ALD_opt ALD_opt];
    ALD_range = [ALD_range ALD_range];
    WTD_opt = [WTD_opt WTD_opt];
    WTD_range = [WTD_range WTD_range];
    ag_frac_npp = [ag_frac_npp ag_frac_npp];
    bg_frac_npp = [bg_frac_npp bg_frac_npp];
    k_0 = [k_0 k_0];
    k_month_0 = [k_month_0 k_month_0];
    NPP_rel = [NPP_rel NPP_rel];
    mosses =    [ mosses mosses ];
    vasculars = [ vasculars vasculars ];
    sedges =    [ sedges sedges ];
    woody =     [ woody woody ];
%     tf_sedge_root = [tf_sedge_root tf_sedge_root];
%     tf_non_sedge_root = [tf_non_sedge_root tf_non_sedge_root];
%     tf_moss = [tf_moss tf_moss];
%     tf_ag_pft = [tf_ag_pft tf_ag_pft];
end


% ****************************
%  Radiocarbon
% ****************************
% READ IN ATMOSPHERIC DEL-14C AND INTERPOLATE TO ANNUAL atm_c14_infile_name_csv = 'atm_del_14C_20,000BP_to_2100_AD_RCP85_RCP45.csv';
atm_c14_infile_name_csv = 'HPM-20 input files\atm_del_14C_20,000BP_to_2500_AD_RCP85_RCP45.csv'; % atm del14-C held constant for 2100-2500 (not realistic)
% atm_c14_infile_name_csv = 'atm_del_14C.csv';
atm_del_c14_time_series = importdata(atm_c14_infile_name_csv);

atm_del_c14_years_input = atm_del_c14_time_series.data(:,1);
atm_del_c14_values_RCP85_input = atm_del_c14_time_series.data(:,2);
atm_del_c14_values_RCP45_input = atm_del_c14_time_series.data(:,3);

% TRUNCATE DEL-C14 TO SIMULATION PERIOD
% sim_start is start year (in years BP; 1950 = 0 BP)
% sim_end is end year (in years BP; 1950 = 0 BP)
% sim_len is simulation length in years (='start' - 'end' + 1 (for 1950 = 0 BP)

index_start = find(atm_del_c14_years_input <= sim_start,1);
index_end = find(atm_del_c14_years_input <= sim_end,1);
atm_del_c14_years = atm_del_c14_years_input(index_start:index_end);
atm_del_c14_values_RCP85 = atm_del_c14_values_RCP85_input(index_start:index_end);
atm_del_c14_values_RCP45 = atm_del_c14_values_RCP45_input(index_start:index_end);

atm_del_c14_simyears = max(atm_del_c14_years) - atm_del_c14_years + 1;  % convert 14C year data to simulation year #

atm_del_c14_ann = interp1(atm_del_c14_simyears, atm_del_c14_values_RCP85, 1:1:sim_len,'pchip');
% atm_del_c14_ann = interp1(atm_del_c14_years, atm_del_c14_values_RCP45, 1:1:sim_len,'pchip');
atm_c14_ann = atm_del_c14_ann / 1000 + 1.;  % convert del-value to 'absolute' value

tau_c14 = 5730. / log(2); % radiocarbon decay rate ('log' in matlab is 'ln')


% ****************************
%  Save parameters
% ****************************

save('hpm20_param_vals','outname', 'climForcingFile', 'sim_len','sim_start','sim_end','tau_c14',...
    'num_veg','ag_frac_npp','bg_frac_npp',...
    'wtd_0','ald_0', 'start_depth', 'lag_years','gipl_flag','gfdl_model_flag','flag_RCP', ...
    'rootin_d80','rootin_alpha','rootin_c3','rootin_c4','rootin_c5',...
    'k_0','k_month_0','wfps_opt','wfps_curve','wfps_sat_rate','wfps_min_rate','anoxia_scale_length',...
    'min_bulk_dens','del_bulk_dens','dens_c1','dens_c2','OM_dens',...
    'ann_temp','ann_temp_orig','ann_temp_amp','ann_temp_offset','ann_ppt', ...
    'ET_0','ET_wtd_1','ET_wtd_2','ET_wtd_3','ET_param','ET_temp_sens',...
    'Roff_c1','Roff_c2','Roff_c2a','Roff_c3','Roff_c4','runon_c1','runon_c2','runon_c3',...
    'wfps_c1','wfps_c2','wfps_c3','mosses','vasculars','sedges','woody',...
    'WTD_opt','WTD_range','PD_opt','PD_range','ALD_opt','ALD_range','NPP_rel','NPP_Q10','max_npp',...
    'tf_old_new','year_old2new','pf_flag');
%     'tf_moss', 'tf_sedge_root', 'tf_non_sedge_root','tf_ag_pft','tf_old_new','year_old2new');
%       'temp_flag','temp_rand_persist','temp_amp1','temp_amp2',...'ppt_flag','ppt_amp1','ppt_amp2','ppt_rand_persist','mon_ppt_offset',...
