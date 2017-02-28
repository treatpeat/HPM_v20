% This version reads in monthly temperature and precipitation from a file;
%   monthly temp used as driver for GIPL2; monthly precip aggregated to annual at this time

% HPM12_soiltemp_gipl2 (winter 2014-15): This version includes old/new carbon switch for arbitrary number of PFTs
%      (requires modifying hpm12_params12.m file)
%   goal is to have single version that works both for permafrost and non-pf
%       - not sure if it is there yet.

% HPM 12 (summer 2014) adding soil temperature routine from GIPL2 model (already implemented in
% HPM12_PF_gipl).

% HPM 10 (latest modifications: April 2012)

% added 14C labeling Oct. 2013

% Model script for Holocene Peatland Model (HPM v.10) (revised PDM model lincked to PAM) 

% HPM v.10 has climate (precipitation, temperature) in a seperate file and tree PFT. 
% HPM v.10 has new routine for root inputs - equally partitioned among all
%    layers in the root zone for non-sedge, and among all layers in the upper
%    (80%) zone for sedges, rather than proportional to layer thickness. This
%    seems to minimize WT oscillations that appeared in runs with constant
%    precipitation and no mosses (Aug. 2011; S. Frolking)
% HPM v.10 has a new algorithm to keep track of moss fraction of peat through simulation.  
%    Max potential peat height 5m (for now), bins of 2cm (for now) (Nov.2011; S. Frolking).
% HPM v.10 does not have WatBal10a routine anymore


% March 2010
% HPM v.9 revised peat ht impact on runoff (more on this still to be done) so that
%    runoff_factor = 1 + params.Roff_c2 * (sum(THICK) - params.Roff_c2a))
% i.e., negative impact until some peat height reached (peats initialize in low areas)
% revised anoxia scale length factor to be proportional to relative
% fraction of NPP due to sedge PFTs  (NOT IMPLEMENTED YET)
% added a new age-weighted mass factor to look at bias in dating due to roots

% October 2009
% HPM v.8 changes to veg NPP and lagWTD based on comments from Jill, Eevastiina, and Lisa B at 
% Prague meeting.  changed MB precip history to more closely match reconstruction in figure of 
% apparent C accumulation. gave all current subroutine version #8.

% HPM v.7 (skipped this one)

% HPM v.6 has changes to water balance subroutine.  It comes much closer to conserving water, 
% by finding a water table depth that is consistent with the new water content (= old water 
% content + net water added); it does not use a specific yield approach.

% HPM v.5 removes 'fast' litter pool (unused in v.4), dropped 'slow' from names. moves all(?) 
% parameters into 'hpm_params5.m'. increases vegetation number and makes a 2-D veg space 
% (Water Table Depth and Peat Depth/ombrotrophy): see file 'vegetation_NPP_2.m'.
% adds WFPS control on decomposition to slow it down under dry conditions.

% HPM v.4 combines all water balance functions into a single function, and adds a transmissivity 
% profile for runoff and a specific yield profile and run-on. uses only the 'slow' mass pool.

% HPM v.3 expands on HPM v.2 by using the PAM formulation to predict productivity and water balance.
   
% HPM v.2 expands on HPM v.1 by including multiple vegetation types (up to 6): woody evergreen, 
% woody deciduous, herb, sedge, sphagnum, non-sphagnum bryophyte.
   
% the revised version (HPM v.1) differs in several respects:
%   1. it builds the peat profile from oldest to youngest (forward time)rather than youngest
%      to oldest (backward time)
%   2. Version 1 will be the basic structure (only one plant type)
%   3. Litter will be partitioned into fast (dm/dt = -k_fast * m) and 
%      slow (dm/dt = -k_slow * (m/m_0) * m
%   4. NPP, WTD, and temp effects will be prescribed (to come from PAM later)
%   5. anoxia transition will be a smooth function (following Bauer 2004)
%   6. root inputs will follow Bauer (2004)

% reference for PDM is:
% Frolking S, NT Roulet, TR Moore, PJH Richard, M Lavoie, & SD Muller. 2001.
%    Modeling northern peatland decomposition and peat accumulation, Ecosystems, 4:479-498.

% reference for PAM is:
% Hilbert D, NT Roulet, TR Moore. 2000. Modelling and analysis of peatlands as dynamical systems, 
% J. Ecology. 88:230-242.

% functions called:
%   hpm_params20
%   hpm_vegNPP20  (called from hpm_params20)
%   hpm_climate_params20 (called from hpm_params20)
%   hpm_npp0
%   hpm_dens0
%   hpm_WatBal0
%   hpm_decomp0
%   hpm_rootin0
%   hpm_climate20  
%   hpm20_gipl2_daily    UAF GI permafrost model from Marchenko
%   hpm_gipl_params20       parameters for UAF GI permafrost model

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% INITIALIZE THINGS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% load in HPM parameters & initialize
% check parameter file, but typical mass units are kg/m^2 dry mass & m water depth (ET,PPT, Runoff, ...)

% hpm_params12_MB;
% hpm_params13x_MB;
hpm_params20x_TLK;
params=load('hpm20_param_vals');

if (params.tf_old_new > 0.5)  % using Old-New carbon switch
    old_new_ones = ones(1,num_veg/2);
    old_new_zeros = zeros(1,num_veg/2);
end


% paramaters for UAF GIPL 2.0 soil temperature model
hpm_gipl_params20;  
params_gipl = load('hpm20_gipl_param_vals');   

nveg = params.num_veg;

num_years = params.sim_len;

% For using Muller MB precip reconstruction, num_years must be multiple of 250 y
% num_years = params.sim_len + 250 - mod(params.sim_len,250)

timestep = 1;   % [y]  BE CAREFUL ABOUT CHANGING THIS FROM ONE (1)!!!
istep = num_years / timestep;

base_ppt = params.ann_ppt;  % annual ppt (m/y) from Roulet PAM

% preallocate arrays to speed up simulations

% small m arrays are masses as annual cohort by veg types
m = zeros(istep,nveg);   % remaining mass in cohort (layer) i and veg type
m_0 = zeros(istep,nveg); % total input mass in cohort i and veg type
m_0_age = zeros(istep,nveg); % total input mass in cohort i and veg type
m_star = zeros(istep,nveg);   % = m / m_0
c14_m = zeros(istep,nveg);  % 14-C in each PFT of each cohort

% capital M vectors are masses as annual cohort accumulated across the veg types
M = zeros(istep,1);        % = sum m across veg. types in cohort/layer i
M_0 = zeros(istep,1);      % = sum m_0 across veg. types in cohort/layer i
M_0_age = zeros(istep,1);      % = sum m_0 across veg. types in cohort/layer i
M_star = zeros(istep,1);        % = M / M_0 in cohort/layer i
M_overlying = zeros(istep,1);   % = sum M_total in profile above cohort/layer i
del_M_tot = zeros(istep,1);   % annual change in total peat mass
c14_M = zeros(istep,1);    % 14-C in each cohort

% these are temporary arrays
mstemp = zeros(istep,nveg);
ms0temp = zeros(istep,nveg);
ms0agetemp = zeros(istep,nveg);
mstartemp = zeros(istep,nveg);
ktemp = zeros(istep,nveg);
agebiastemp = zeros(istep,nveg);
agebiastemparr = zeros(istep,nveg);

% vectors down the profile
depth = zeros(istep,1);  % cohort (layer) depth in meters
thick = zeros(istep,1);  % cohort (layer) thickness in meters
zbottom = zeros(istep,1); % depth (m) from top of peat to bottom of cohort
porosity = zeros(istep,1);  % cohort (layer) porosity (m3/m3)
prev_thick = zeros(istep,1);  % cohort (layer) thickness in meters (from previous time step)
dens = zeros(istep,1);   % cohort (layer) bulk density in kg/m3
time = zeros(istep,1);   % keps track of time in years
age_bias = zeros(istep,1);  % for keeping track of age bias
tmp_depth = zeros(500,1); % temporary truncated array
profiles = zeros(istep,6*floor(params.sim_len - params.year_old2new)/50.);  % write 6 variable profile every 50 years near end of simulation

% these are temporary arrays
depth2 = depth;
dens1 = dens;
dens_old = dens;
dens_old2 = dens_old;
dens_evolve = zeros(istep,4);

% arrays by cohort and veg type
k = zeros(istep,nveg);  % mass loss rate (1/y)

% vectors down the profile
k_mean = zeros(istep,1);        % mass-weighted mean decomposition factor by cohort
anoxiafact = zeros(istep,1);    % anoxia profile, function of water table depth (anything else?)

% array of root mass input (kg/m2/layer) by veg type
rootin = zeros(istep,nveg);
rootin2 = zeros(istep,nveg);  % temporary array

% arrays by time and veg type
annNPP = zeros(istep,nveg);
tot_npp = zeros(istep,1);
annRESP = zeros(istep,1);  % annual mass loss (carbon units = biomass/2)
annROOTIN = zeros(istep,1);
annROOTNPP = zeros(istep,1);
annAGMASSIN = zeros(istep,1);
annZ_total = zeros(istep,1);
del_peat_height = zeros(istep,1);
annM_total = zeros(istep,1);
del_c14_annRESP = zeros(istep,1); % del_14c of annual respiration

% NPPVEC = zeros(nveg);
del_C_del_t = zeros(istep,1);
del_C_del_t2 = zeros(istep,1);
j5 = zeros(istep,1);
atm_del_14C = zeros(istep,1);  % time series of atmos. del-14C (read in from file)

% vectors and arrays for debugging, etc.
junk1 = zeros(istep,3);
junk2 = zeros(istep,3);
temporary = zeros(istep,1);
cohortM = zeros(istep,10);

% vectors by time 
annPPT = zeros(istep,1);   % m/yr
annWTD = zeros(istep,1);   % m
peat_water = zeros(istep,1);  % m3/m2
total_water = zeros(istep,1);  % m3/m2
lagWTD = zeros(istep,1);  % years
annTRANS = zeros(istep,1);   % relative hydraulic transmissivity (0-1)
WATER = zeros(istep,7);   % array for output that contains annual water balance terms
annTHETA = zeros(istep,1);
annWTD_VAR = zeros(istep,1);
annWFPS = zeros(istep,1);
% prev_annWFPS = zeros(istep,1);
del_peatwater = zeros(istep,1);
net_water_in = zeros(istep,1);
annTEMP_FACT = ones(istep,1);
ann_ET_Tfact = zeros(istep,1); % temperature impact on annual ET (multiplier)
counter_array = zeros(istep*10,2);
reconstrWTD = zeros(istep,1);
age_depth = -9999 * ones(istep,30);
age_depth2 = -9999 * ones(istep,375);
annTEMP = ones(istep,1);
annTEMP_NPP_FACT = ones(istep,1);  % NPP temperature modifier

% GIPL SOIL TEMP OUTPUT
ALD1_gipl = zeros(istep,1);  % GIPL model annual max active layer depth to Tfr + FIT
ALD2_gipl = zeros(istep,1);  % GIPL model annual max active layer depth to Tfr
ALD3_gipl = zeros(istep,1);  % GIPL model annual max active layer depth to Tfr - FIT
soil_node_temp_month_save = zeros(istep,12,params_gipl.NumberOfSoilComputationNodes);
Ann_snow_depth = zeros(istep,1);  % GIPL model max annual snowdepth (meters?)
ann_ALD = zeros(istep,1);  % annual active layer max depth (Kudryavtsev)

% vectors and arrays for the math
onevec = ones(istep,1);
epsvec = eps*ones(istep,1);
zerovec = zeros(istep,1);
onearr = ones(istep,nveg);
epsarr = eps*ones(istep,nveg);
topvec = zeros(1,nveg);
topval = 0;

% post-disturbance old and new carbon variables   % modified for HPM10PF
M_old = zeros(istep,1);        % = sum m across 'old' veg. types in cohort/layer i
M_new = zeros(istep,1);        % = sum m across 'new' veg. types in cohort/layer i
ann_M_old = zeros(istep,1);    % annual total 'old' M
ann_M_new = zeros(istep,1);    % annual total 'new' M
ann_NPP_old = zeros(istep,1);    % annual total 'old' NPP
ann_NPP_new = zeros(istep,1);    % annual total 'new' NPP
cohort_age = zeros(istep,1);           % root-input-adjusted cohort age
cohort_age_temp = zeros(istep,1);           % root-input-adjusted cohort age
ann_resp_age = zeros(istep,1);   % age-mass weight of annual respiration
cohort_age2 = zeros(istep,nveg);           % root-input-adjusted cohort age
cohort_age_temp2 = zeros(istep,nveg);           % root-input-adjusted cohort age
ann_resp_age2 = zeros(istep,1);   % age-mass weight of annual respiration
ann_resp_old = zeros(istep,1);    % annual total 'old' resp/decomp
ann_resp_new = zeros(istep,1);    % annual total 'new' resp/decomp

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% END OF INITIALIZATIONS
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% *********** moss fraction *********************************************
% (S. Frolking)
% variables for binning M* and moss fraction of peat ********

flag_bins = 1;   %  If flag_bins = 1, compute M*, if 2, compute M* & moos_frac

if (flag_bins > 0.5) 
    nbins = 250;                % for binning cohorts in output
    maxheight = 3.; % 6.;             % total potential height (meters) 
    delx = maxheight/nbins;     % total possible ht (meters) ÷ # of bins
    cohortheight = zeros(istep,1);      % height of top of cohort above bottom of peat
    bin_M_star = 9999 * ones(nbins,istep); % bin M* = M/M_0
end

if (flag_bins > 1.5) 
    mossfrac = zeros(istep,1);          % cohort mass fraction that is moss
    bin_moss_frac = -0.9999 * ones(nbins,istep); % bin moss fraction that is moss
end

% *****************************************************

% initialize surface cohort with aboveground litter inputs from all plant types

time(1) = timestep / 2 ;
thick(1) = 0.05;  % placeholder value for first year NPP calculation

NPP = hpm_npp20(params.wtd_0, params.wtd_0, params.ald_0, thick, params);  % use permafrost version

if (params.tf_old_new > 0.5)
    NPP = NPP  .* [old_new_ones old_new_zeros];   % modified for 2xN PFTs
end

m(1,:) = NPP .* params.ag_frac_npp;
m_0 = m;
m_0_age = m * (num_years - 0.5);
m_star = m ./ (epsarr + m_0);
age_bias(1) = 1;
c14_m(1,:) = m(1,:) * atm_c14_ann(1);

M = sum(m,2);
M_0 = M;
M_0 = M * (num_years - 0.5);
M_star = M ./ (epsvec + M_0);
c14_M(1) = sum(m(1,:) .* c14_m(1,:)) / sum(m(1,:));

M_overlying(1) = 0;
prev_M_tot = 0;

annNPP(1,:) = NPP;
annWTD(1) = params.wtd_0;
annWTD_VAR(1) = annWTD(1)/3;
annTEMP_FACT(1) = 1.0;

ann_ALD(1) = 0.25;  % arbitrary first year ALD (m)

% calculate layer density, thickness, and depth

dens = hpm_dens20(M_star, M_overlying, params, onevec);

thick(1) = M(1) / (eps + dens(1));
zbottom(1) = thick(1);
prev_thick(1) = thick(1);
depth = cumsum(thick) - onevec * thick(1)/2;
annZ_total(1) = thick(1);
annM_total(1) = M(1);
% last_Z_total = 0;

tic;

flag1 = 0;   % set to 1 when simulation of dynamic water balance begins
flag2 = 0;   % set to 1 when simulation of dynamic water balance begins
profile_counter = 1;  % used to write out some profiles near end of run

%-- RECORD ANNUAL PRECIP AND TEMPERATURE VALUES 
REC_annppt = precip_forcing_ann;
REC_anntemp = temp_forcing_ann;

mean_REC_annppt = mean(REC_annppt);
REC_annppt_orig = REC_annppt;

%-- STORE SOME PARAMETER VALUES FOR ADJUSTMENT DURING SENSITIVITY RUNS
Roff_c2_orig = params.Roff_c2;
ET_0_orig = params.ET_0;
NPP_rel_orig = params.NPP_rel;
NPP_rel_orig_moss_frac = sum(NPP_rel_orig .* params.mosses)/(eps + sum(NPP_rel_orig));
NPP_rel_vasc_orig = NPP_rel_orig .* params.vasculars;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% LOOP THROUGH YEARS OF SIMULATION
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for itime = 2:timestep:num_years

    time(itime) = (itime - 0.5) * timestep;

    T_mean = REC_anntemp(itime);
    REC_anntemp_amp(itime) = params.ann_temp_amp;

    if (mod(itime,500*timestep) == 0)   % tracks/writes out clock time per 500 y of simulation
        toc;
        tic;
        timex = (itime - 0.5) * timestep
    end
    
    if (itime < 2.5)
        delpeat = annZ_total(itime-1);
    else
        delpeat = annZ_total(itime-1) - annZ_total(itime-2);
    end
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% CALL DRYING/DRAINAGE, WARMING OR OTHER CLIMATE SCENARIOS HERE
% load('hpm_dryingScenarios.m')
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% CALCULATE ANNUAL WATER BALANCE IN MONTHLY TIME STEPS 

% net water increase = precipitation + runon - (evapotranspiration + runoff)
    
    if (depth(itime-1) < params.start_depth)     % don't use dynamic water table until there is enough peat
        initflag = 1;
        annWTD(itime) = params.wtd_0; 
        lagWTD(itime) = params.wtd_0;
        annWFPS = 0.8 * onevec;
        WATER(itime,:) = [0 annPPT(itime) 0 0 0 annWTD(itime) delpeat*100];
        junk1(itime,:) = [0 0 0];
%        junk2(itime,:) = [anndelwat annPPT(itime) annRUNON-annET-annRUNOFF];
    
    else
        initflag = 0;
        if (flag1 <0.5)
            flag1 = 1   % only do this stuff once  ??? or do each time 'leaving' depth < 0.15 m condition???

            monWTD = params.wtd_0;
            annWTD(itime-1) = monWTD;
            
            if (monWTD >= max(depth))   % if WT is below peat, set index to near bottom of peat
                index = time - 1;
            else
                index = find(depth>monWTD,1);   % finds vector index value of first case of layer depth > WTD
            end

            porosity = onevec - dens/params.OM_dens;  % bulk density = mass peat ÷ (volume organic matter + porosity)

            zstar = params.wfps_c1 * onevec + (params.wfps_c2 - params.wfps_c1)*((dens - params.min_bulk_dens)...
                    ./(dens - params.min_bulk_dens + params.wfps_c3));
     
            zwtd = depth - monWTD;  % determines distance each cohort is from WT (value is positive if cohort is below WT, i.e., submerged)
            zwtd = max(zerovec, -zwtd);    % determines distance above WT, set to zero if at or below WT

            annWFPS = params.wfps_c1 + (1 - params.wfps_c1) * exp(-zwtd./zstar);   % see notes and file 'anoxia & bulk dens & WFPS % profile.xls')

            new_peat_wat = sum((annWFPS .* thick) .* porosity);
            flag1reportpeatwater = new_peat_wat
            
            dynamic_watbal_time_start = itime
            peat_depth = annZ_total(itime-1)
            del_ppt = 0;  % ???
            cohortmass = M(1:itime)
        end
        
       % annual precipitation
        
      annPPT(itime) = REC_annppt(itime) + annppt_modifier;

       %%%%%%%%%%%%%%%%%%%%%%%%%%     
      
        
        stepperyr = 3;  % # subannual water balance time steps (could become seasonal water balance)
        
        if (flag2 < 0.5)
            flag2 = 1;
            monRUNON = 0;  % for first dynamic water balance step
            monRUNOFF = 0; % for first dynamic water balance step
            monET = 0.7 * (params.ann_ppt / stepperyr);   % for first dynamic water balance step 
        end            
 
% zero annual water balance values, then compute multiple steps per year to get mean value
        anndelwat = 0;
        annET = 0.;
        annRUNOFF = 0.;
        annRUNON = 0.;
        annTRANS(itime) = 0;
        ann_wtd = 0;
        annWFPS = zerovec;

        monWTD = annWTD(itime-1);
        total_porosity = sum(thick .* porosity);   % total peat porosity = water content at saturation

        zstar = params.wfps_c1 * onevec + (params.wfps_c2 - params.wfps_c1)*((dens - params.min_bulk_dens)...
                 ./(dens - params.min_bulk_dens + params.wfps_c3)); % water content as a function of decomposition
             
        zbottom = cumsum(thick);
        
% zstar is a parameter for WFPS calculation, and is a function of bulk density (i.e., m/m0)
     
        old_peat_wat = new_peat_wat;
        
        for ii = 1:1:stepperyr
            monPPT = annPPT(itime) / stepperyr;
            delwat = monPPT + monRUNON - monRUNOFF - monET;
            anndelwat = anndelwat + delwat;
            peat_wat = new_peat_wat;
            total_wat = peat_wat + max(0, -monWTD);
            del_peat_ht = 0;
            
% new water balance calculation (June 2009)
%   INPUTS: WTD, del_water, peat_water, density, thick, porosity, total_poros, zstar, ones, zeros, params, time
%  OUTPUTS: new_WTD, new_water_content, new_peat_water, ET, R_off, R_on, transmissivity

            if (abs(delwat) > 0.0000001)
%                 [new_wat_cont, new_wtd, new_peat_wat, annet, montrans,annrunoff,annrunon, counter1, counter2] = ...
%                     hpm_WatBal10(T_mean, del_peat_ht, monWTD, delwat, peat_wat, dens, thick, depth, porosity, total_porosity, zstar, zbottom,onevec, zerovec, params,(itime-1));

                if params.pf_flag < 0.5  % NO PERMAFROST
                    [new_wat_cont, new_wtd, new_peat_wat, annet, montrans,annrunoff,annrunon, counter1, counter2] = ...
                        hpm_WatBal20_test(T_mean, del_peat_ht, monWTD, delwat, peat_wat, dens, thick, depth, porosity, ...
                            total_porosity, zstar, zbottom,onevec, zerovec, params,(itime-1));
                else
                    [new_wat_cont, new_wtd, new_peat_wat, annet, annet_temp_fact, montrans,annrunoff,annrunon,deep_count] = ...
                        hpm_WatBal20pf(monWTD, delwat, peat_wat, dens, thick, depth, porosity, total_porosity, zstar, zbottom, ...
                            onevec,zerovec, params,(itime-1),ann_ALD(itime-1),T_mean);
                end
            end
   
%             if (abs(counter1 - counter2) > 100)
%                 search_steps = abs(counter1-counter2)
%                 year = itime
%                 month = ii
%             end
    
            new_tot_wat = new_peat_wat + max(0, -new_wtd);
            
            if ((new_peat_wat - (peat_wat + delwat)) > delwat)
                result_NewPeatWat_OldPeatWat_delwat_year_timestep = [new_peat_wat peat_wat delwat*100 itime ii];
                result_monppt_monrunon_monrunoff_monET_newwtd = [monPPT monRUNON monRUNOFF monET new_wtd];
            end
  
            monET = annet / stepperyr;       % units of flows in hpm_WatBal6 are m/y; this makes m/timestep
            monRUNOFF = annrunoff / stepperyr;
            monRUNON = annrunon / stepperyr;

            annET = annET + monET;
            annRUNOFF = annRUNOFF + monRUNOFF;
            annRUNON = annRUNON + monRUNON;
            annTRANS(itime) = annTRANS(itime) + montrans / stepperyr;
            annWFPS = annWFPS + new_wat_cont / stepperyr;

            monWTD = new_wtd;
            ann_wtd = ann_wtd + monWTD / stepperyr;
            
        end
        annWTD(itime) = ann_wtd;   % mean value for 'stepperyr' sub-annual time steps
%        ann_ET_Tfact(itime) = annET_Tfact;
        del_peatwater(itime) = new_peat_wat - old_peat_wat; % annual change in peat water content
        net_water_in(itime) = anndelwat;   % annual net water input
        peat_water(itime) = new_peat_wat;  % end of year value
        total_water(itime) = new_tot_wat;  % end of year value
        
    end   % END OF ANNUAL WATER BALANCE CALCULATION 

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% CALCUALTE TEMPERATURE EFFECTS AND ACTIVE LAYER DEPTH
%    using daily timestep GIPL2 model (hpm12_gipl2_daily.m)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
%  --- start by generating a monthly time series from the annual temp and annual temp amplitude

    if(itime == 2)

        T_jan1_gipl = params_gipl.T_init_gipl;
%        T_jan1_gipl = params_gipl.T_init_gipl - 10.;
        
        swe_prev = 0.1;
        
        months = 13:1:24;  % note that model starts with year 1 cohort in place
    else 
        months = months + 12;
    end
    
    air_temp_month = temp_forcing_mon(months);
    monthly_precip = precip_forcing_mon(months);  
    
%  --- use GIPL2 model at daily timestep to compute:
%       monthly GIPL soil layer temps ('params_gipl.ndepth' layers to 10m, excludes bedrock)
%       monthly GIPL soil node temps (to 100 m)
%       final year day snow-water equiv (for initializing next year)
%       final year day soil temp profile (for initializing next year)
%       annual max GIPL active layer (for Tfr+FIT; Tfr; and Tfr-FIT)

    [soil_layer_temp_month soil_node_temp_month swe_day_prev soil_temp ALD1 ALD2 ALD3 snow_depth_max] = ... 
        hpm20_gipl2_daily(itime, T_jan1_gipl, air_temp_month, monthly_precip, swe_prev, depth, thick, annWFPS, porosity, dens, params, params_gipl);
    
% save Dec profile to initiate next year
    swe_prev = swe_day_prev;
    T_jan1_gipl = soil_temp;
%    T_init_gipl = soil_layer_temp_month(12,:);   % save Dec profile to initiate next year

    ALD1_gipl(itime) = ALD1;
    ALD2_gipl(itime) = ALD2;
    ALD3_gipl(itime) = ALD3;
    
   if params.pf_flag > 0.5
       ann_ALD(itime) = max(0.15,ALD2_gipl(itime));   % if ALD = 0, root input becomes NaN, then mass becomes NaN, then crash
   else
       ann_ALD(itime) = 200. % put active layer below everything
   end

    soil_node_temp_month_save(itime,:,:) = soil_node_temp_month;
    Ann_snow_depth(itime) = snow_depth_max;
    annTEMP(itime) = mean(air_temp_month);

%    annTEMP_NPP_FACT(itime) = 1.5 ^ ((T_mean - params.ann_temp_orig)/10); % Q10 = 1.5, ref temp = base temp.    
    annTEMP_NPP_FACT(itime) = params.NPP_Q10 ^ ((T_mean - params.ann_temp_orig)/10); % Q10 = 1/1.5 = 0.67, ref temp = base temp.
 
% interpolate monthly GIPL2 node temperatures to HPM cohort depths

    T_month_cohort = interp1(soilLayerDepth',soil_node_temp_month', depth');
    T_month_cohort = T_month_cohort';
    
%  WTD lag time for NPP   
    
    if (itime < params.lag_years + 1)
        lagWTD(itime) = mean(annWTD(1:itime));
    else
        lagWTD(itime) = mean(annWTD(itime-params.lag_years:itime));
    end
    
% calculate annual productivity for each veg type as function of WTD

    NPP = hpm_npp20(annWTD(itime),lagWTD(itime), ann_ALD(itime), thick, params);

    NPP = NPP * annTEMP_NPP_FACT(itime);

    if (params.tf_old_new > 0.5) % if tf_old_new = 1, then implement old/new PFTs
           
        if (itime < params.year_old2new)
            NPP = NPP .* [old_new_ones old_new_zeros];  % use 'old' PFTs 
        else
           NPP = NPP .* [old_new_zeros old_new_ones];  % use 'new' PFTs 
        end
        
    end

    annNPP(itime,:) = NPP(:);
    tot_npp(itime) = sum(NPP);
    if (params.tf_old_new > 0.5) 
        ann_NPP_old(itime) = sum(NPP(1:num_veg/2));
        ann_NPP_new(itime) = sum(NPP(num_veg/2+1:num_veg));
    end

% determine root inputs
        
%    rootin = hpm_rootin10(depth, thick, params, NPP, annWTD(itime), annZ_total(itime-1), onevec);
    rootin = hpm_rootin20(depth, thick, params, NPP, annWTD(itime), ann_ALD(itime-1), annZ_total(itime-1), onevec);
             % hpm_rootin12 is for arbitray number of PFTs
             
    rootin2 = [topvec; rootin];
    rootin2(end,:) = [];
    rootin = rootin2;
    annROOTIN(itime) = sum(sum(rootin,2));
    annROOTNPP(itime) = sum(NPP .*params.bg_frac_npp);
    j5(itime) = annROOTIN(itime) - sum(NPP .* params.bg_frac_npp); %debug

% store masses in temporary vectors, adding zero to top, then removing final (bottom) zero to maintain same total size
%   repeat this routine a few lines later decomp array, k, and on the root input array

    mstemp = [topvec; m];    % add zeros to top row
    ms0temp = [topvec; m_0];
    ms0agetemp = [topvec; m_0_age];
    mstartemp = [topvec; m_star];
    agebiastemp = [topval; age_bias];
    c14_mstemp = [topvec; c14_m];    % add zeros to top row

    mstemp(end,:) = [];         % remove final row (of zeros) to maintain array size
    ms0temp(end,:) = [];
    ms0agetemp(end,:) = [];
    mstartemp(end,:) = [];
    agebiastemp(end) = [];
    agebiastemparr = repmat(agebiastemp,1,nveg);
    c14_mstemp(end,:) = [];
    
    depth2 = [topval; depth];
    depth2(end) = [];
    depth = depth2;

% determine factors to age cohorts by one year

% DETERMINE DECOMPOSITION

% -- DECOMPOSITION ANNUAL WATER/ANOXIA MULTIPLIER

% 'decompfact' returns annual profile of decomp rate modifier for peat water content and WT position
    decompfact_water = hpm_decomp20(depth, annWTD(itime), annWFPS, params, onevec, epsvec);
    
% -- DECOMPOSITION MONTHLY TEMPERATURE MULTIPLIER
% temp effect is continuously varying Q10, with 10�C reference 
%    see notes in 'annual vs. monthly time step decomposition time step 1.xlsx'
%    Q10(T) = 2 + 3 * exp(-T/9), with T in celsius degrees 
%     this follows Q10 curve from Hamdi et al. 2013, and is similar to Lloyd & Taylor function used
%     by Wania et al. 2009 & Spahni et al. 2013.
%    params_gipl.Tfr = midpoint of phase change range (degC); params.gipl.FIT = half-width of phase change range (degC)

    if (params.gipl_flag > 0.5)
%         decompfact_temp = (T_month_cohort > (params_gipl.FIT + params_gipl.Tfr)) .* ...
        decompfact_temp = (T_month_cohort > (params_gipl.Tfr)) .* ...
            (2 + 3*exp(-T_month_cohort/9)).^((T_month_cohort - 10)/10);
        decompfact_temp = decompfact_temp';
        k_jan = (decompfact_water .* decompfact_temp(:,1) * params.k_month_0) .* mstartemp;
        k_feb = (decompfact_water .* decompfact_temp(:,2) * params.k_month_0) .* mstartemp;
        k_mar = (decompfact_water .* decompfact_temp(:,3) * params.k_month_0) .* mstartemp;
        k_apr = (decompfact_water .* decompfact_temp(:,4) * params.k_month_0) .* mstartemp;
        k_may = (decompfact_water .* decompfact_temp(:,5) * params.k_month_0) .* mstartemp;
        k_jun = (decompfact_water .* decompfact_temp(:,6) * params.k_month_0) .* mstartemp;
        k_jul = (decompfact_water .* decompfact_temp(:,7) * params.k_month_0) .* mstartemp;
        k_aug = (decompfact_water .* decompfact_temp(:,8) * params.k_month_0) .* mstartemp;
        k_sep = (decompfact_water .* decompfact_temp(:,9) * params.k_month_0) .* mstartemp;
        k_oct = (decompfact_water .* decompfact_temp(:,10) * params.k_month_0) .* mstartemp;
        k_nov = (decompfact_water .* decompfact_temp(:,11) * params.k_month_0) .* mstartemp;
        k_dec = (decompfact_water .* decompfact_temp(:,12) * params.k_month_0) .* mstartemp;
    else
        decompfact_temp = ones(istep,12);
        k_jan = (decompfact_water .* decompfact_temp(:,1) * params.k_0 / 12) .* mstartemp;
        k_feb = (decompfact_water .* decompfact_temp(:,2) * params.k_0 / 12) .* mstartemp;
        k_mar = (decompfact_water .* decompfact_temp(:,3) * params.k_0 / 12) .* mstartemp;
        k_apr = (decompfact_water .* decompfact_temp(:,4) * params.k_0 / 12) .* mstartemp;
        k_may = (decompfact_water .* decompfact_temp(:,5) * params.k_0 / 12) .* mstartemp;
        k_jun = (decompfact_water .* decompfact_temp(:,6) * params.k_0 / 12) .* mstartemp;
        k_jul = (decompfact_water .* decompfact_temp(:,7) * params.k_0 / 12) .* mstartemp;
        k_aug = (decompfact_water .* decompfact_temp(:,8) * params.k_0 / 12) .* mstartemp;
        k_sep = (decompfact_water .* decompfact_temp(:,9) * params.k_0 / 12) .* mstartemp;
        k_oct = (decompfact_water .* decompfact_temp(:,10) * params.k_0 / 12) .* mstartemp;
        k_nov = (decompfact_water .* decompfact_temp(:,11) * params.k_0 / 12) .* mstartemp;
        k_dec = (decompfact_water .* decompfact_temp(:,12) * params.k_0 / 12) .* mstartemp;
    end

    m_p1 = mstemp .* (onearr - k_jan);
    m_p2 = m_p1 .* (onearr - k_feb);
    m_p3 = m_p2 .* (onearr - k_mar);
    m_p4 = m_p3 .* (onearr - k_apr);
    m_p5 = m_p4 .* (onearr - k_may);
    m_p6 = m_p5 .* (onearr - k_jun);
    m_p7 = m_p6 .* (onearr - k_jul);
    m_p8 = m_p7 .* (onearr - k_aug);
    m_p9 = m_p8 .* (onearr - k_sep);
    m_p10 = m_p9 .* (onearr - k_oct);
    m_p11 = m_p10 .* (onearr - k_nov);
    m_p12 = m_p11 .* (onearr - k_dec);
    
    m = m_p12;
    c14_m = c14_mstemp .* m_p12 ./ (mstemp + epsarr);
    m = m + rootin;
    c14_m = c14_m + rootin * atm_c14_ann(itime);  % remove decomp. and add root c14
    c14_m = c14_m * exp(-1/params.tau_c14);    % annual loss to radioactive decay
        
    ann_resp_array = (mstemp - m_p12) * 0.5 ;  % 0.5 factor for biomass to C
    
    if (params.tf_old_new > 0.5)
        ann_resp_new(itime) = sum(sum(ann_resp_array(:,num_veg/2+1:num_veg),2));
        ann_resp_old(itime) = sum(sum(ann_resp_array(:,1:num_veg/2),2));
    end
       
    annRESP(itime) = (sum(sum(ann_resp_array)));  
    del_c14_annRESP(itime) = sum(sum(((c14_mstemp ./ (eps + mstemp) - 1) * 1000) .* (mstemp - m_p12),2))/sum(sum((mstemp - m_p12),2));

    m_0 = ms0temp + rootin;
    m_0_age = ms0agetemp + rootin * (num_years - time(itime));
%     age_bias = sum((mstemp .* (onearr - k) .* (agebiastemparr+1) + rootin) ./ (mstemp .* (onearr - k) + rootin + epsarr),2) / nveg;
        
%  add new top cohort
    
    m(1,:) = NPP .* params.ag_frac_npp;
    annAGMASSIN(itime) = sum(annNPP(itime,:) .* params.ag_frac_npp);
    m_0(1,:) = m(1,:);
    m_0_age(1,:) = m(1,:) * (num_years - time(itime));
    age_bias(1) = 1;
    c14_m(1,:) = m(1,:) * atm_c14_ann(itime);
    
%    calculate new cohort density and thickness
    
    m_star = m ./ (epsarr + m_0);

    del_C_del_t(itime) = (annAGMASSIN(itime) + annROOTIN(itime))/2 - annRESP(itime);

    junk2(itime,:) = [sum(annNPP(itime,:)) annROOTIN(itime) sum(annNPP(itime,:).*params.ag_frac_npp)];
    
    M = sum(m,2);
    M_tot = sum(M);
    
    if (params.tf_old_new > 0.5) 
        M_old = sum(m(:,1:num_veg/2),2);    % modified for old/new carbon
        ann_M_old(itime) = sum(M_old);
        M_new = sum(m(:,num_veg/2+1:num_veg),2);    % modified for  old/new carbon
        ann_M_new(itime) = sum(M_new);
    end
    
    del_M_tot(itime) = M_tot - prev_M_tot;
    prev_M_tot = M_tot;
    
    M_0 = sum(m_0,2);
    M_0_age = sum(m_0_age,2);
    M_star = M ./ (epsvec + M_0);
    M_overlying = cumsum(M) - M;

%  calculate cohort densities, thicknesses, and depths

    dens1 = hpm_dens20(M_star, M_overlying, params, onevec);
    dens_old2 = [topval; dens_old];
    dens_old2(end) = [];
    dens_old = dens_old2;
    dens = max(dens1,dens_old);   % don't let bulk density decline (??)
    dens_old = dens;
    porosity = onevec - dens/params.OM_dens;  % bulk density = mass peat ÷ (volume organic matter + porosity)

%    dens_evolve(itime,1) = dens(200);   %to seehow densities of particular cohorts change with time
%    dens_evolve(itime,2) = dens(500);
%    dens_evolve(itime,3) = dens(1000);
%    dens_evolve(itime,4) = dens(2000);
    
    prev_thick = thick;
    thick = M ./ (epsvec + dens);
    zbottom = cumsum(thick);
 
    depth = cumsum(thick) - thick/2;
    total_porosity = sum(thick .* porosity);   % total peat porosity = water content at saturation
    
    annZ_total(itime) = depth(itime)+thick(itime)/2;    %    annZ_total(itime) = sum(thick);   % = depth(itime)+thick(1)/2;
    del_peat_height(itime) = annZ_total(itime) - annZ_total(itime-1);
    annM_total(itime) = sum(M);

  
% *********** moss fraction *********************************************
    % calculate moss fraction of peat in 'nbins' vertical bins over 'maxheight' meters from base 
    %       (can be greater than total peat height; missing value is
    %       -0.9999) (S. Frolking)
    if (flag_bins > 0.5)
        
        cohortheight = flipud(cumsum(flipud(thick))); 
        if (flag_bins > 1.5)
            mossfrac = (m .* (onevec * params.mosses)) ./ (M + eps);
        end

    
        x1 = 0.;
        for ix = 1:1:nbins
        
            if (x1 > max(cohortheight))
                break;
            end
            x2 = ix * delx;
    
            tf_bin = (cohortheight > x1) & (cohortheight <= x2);
            tf_bin_sum = sum(tf_bin);

            if(tf_bin_sum>0)
                bin_M_star(ix,itime) = sum(M_star .* M .* tf_bin) / (sum(M .* tf_bin) + eps);
                if (flag_bins > 1.5)
                    bin_moss_frac(ix,itime) = sum(mossfrac .* M .* tf_bin) / (sum(M .* tf_bin) + eps);
                end
            end
    
            x1 = x2;
        
        end     
          
    end
% *********************************************************************    
    
% Reposition water table given peat & total water content and new peat height and total porosity

    if (initflag < 0.5)  % don't during initialization (i.e., initflag = 1)

        if (abs(del_peat_height(itime)) > 0.000001)  % don't bother for small changes in peat height (<0.001 mm)

            peat_wat = sum((annWFPS .* thick) .* porosity);
            zstar = params.wfps_c1 * onevec + (params.wfps_c2 - params.wfps_c1)*((dens - params.min_bulk_dens)...
                     ./(dens - params.min_bulk_dens + params.wfps_c3));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% replaced call to hpm_WatBal10a with call to hpm_WatBal10; July 2011; SF

%            [new_wat_cont, new_wtd, new_peat_wat] = ...
%              hpm_WatBal10a(annWTD(itime),del_peat_height(itime),peat_water(itime),total_water(itime),dens,thick,depth,porosity,total_porosity,zstar,zbottom,onevec,zerovec,params,itime);
            delwat_zero = 0;
%             [new_wat_cont, new_wtd, new_peat_wat, annet, montrans,annrunoff,annrunon,counter1,counter2] = ...
%                 hpm_WatBal10(T_mean, del_peat_height(itime), annWTD(itime), delwat_zero, peat_water(itime), dens, thick, depth, porosity, total_porosity, zstar, zbottom,onevec, zerovec, params,(itime-1));

            if params.pf_flag < 0.5  % NO PERMAFROST
                [new_wat_cont, new_wtd, new_peat_wat, annet, montrans,annrunoff,annrunon,counter1,counter2] = ...
                hpm_WatBal20_test(T_mean, del_peat_height(itime), annWTD(itime), delwat_zero, peat_water(itime), ...
                    dens, thick, depth, porosity, total_porosity, zstar, zbottom,onevec, zerovec, params,(itime-1));
            else
                hpm_WatBal20apf(annWTD(itime),del_peat_height(itime),peat_water(itime),total_water(itime),...
                    dens,thick,depth,porosity,total_porosity,zstar,zbottom,onevec,zerovec,params,itime,ann_ALD(itime));
            end
            
            if params.pf_flag < 0.5
                if (abs(counter1 - counter2) > 100)
                    search_steps = abs(counter1-counter2)
                    year = itime
                end
            end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
         
            annWFPS = new_wat_cont;  
            annWTD(itime) = new_wtd;      % end of year value
            peat_water(itime) = new_peat_wat;  % end of year value

%            junk1(itime,:) = [anndelwat annPPT(itime)+annRUNON-annET-annRUNOFF 0];
%            junk2(itime,:) = [anndelwat annPPT(itime) annRUNON-annET-annRUNOFF];
%            cohortmass = M(1:itime)
        end

        WATER(itime,:) = [anndelwat*100 annPPT(itime) annET annRUNOFF annRUNON annWTD(itime) delpeat*100];
    end

    for i_age_depth = 1:1:30
        if itime > i_age_depth * 500
            age_depth(itime,i_age_depth) = annZ_total(itime) - depth(itime - i_age_depth*500);
        end
    end
  
    for i_age_depth2 = 1:1:375
        if itime > i_age_depth2 * 50
            age_depth2(itime,i_age_depth2) = annZ_total(itime) - depth(itime - i_age_depth2*50);
        end
    end

    if (mod(itime,50*timestep) == 0 && itime >= year_old2new)   % tracks/writes out clock time per 500 y of simulation
        col1 = 1 + (profile_counter-1)*6;
        col2 = profile_counter*6;
        c14_M = sum(c14_m,2) ./ (M + epsvec);
        del_c14_M = (c14_M - 1) * 1000;   % final del-14C profile

        profiles(:,col1:col2) = [ itime*onevec depth  M  M_old  M_new  del_c14_M ];
                       
        profile_counter = profile_counter + 1;
        
    end
    
end % loop through years

monthly_node_temps = zeros(sim_len*12,114);

for jyear = 1:1:sim_len
    for jmonth = 1:1:12
        monthly_temps((jyear-1)*12 + jmonth, :) = soil_node_temp_month_save(jyear,jmonth,:);
    end
end

monthly_layer_temps_300yr = zeros(300*12,63);
for jyear = 1:1:300
    for jmonth = 1:1:12
        monthly_layer_temps_300yr((jyear-1)*12 + jmonth, :) = 0.5 * ...
            (soil_node_temp_month_save(sim_len-300+jyear,jmonth,1:63)+ soil_node_temp_month_save(sim_len-300+jyear,jmonth,2:64));
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% CALCULATE SOME FINAL METRICS AND WRITE OUT & PLOT RESULTS
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

age = time;
years_BP = params.sim_start - time;
M_TOTAL = sum(M);
M_TOTAL2 = sum(del_C_del_t);
Z_TOTAL = depth(end);

k_mean = sum(m .* k,2) ./ (M + epsvec);
c14_M = sum(c14_m,2) ./ (M + epsvec);
del_c14_M = (c14_M - 1) * 1000;   % final del-14C profile

% reconstWTD = zeros(istep,1);
% reconstWTD(:,1) = (m(:,5) * params.WTD_opt(5) + m(:,6) * params.WTD_opt(6) + m(:,7) * params.WTD_opt(7)...
%     + m(:,8) * params.WTD_opt(8) + m(:,9) * params.WTD_opt(9))...
%     ./ (m(:,5) + m(:,6) + m(:,7) + m(:,8) + m(:,9) + eps);
% 

log10junk = 2.14287*onevec - 0.042857 * dens;
hydrconjunk = exp(log(10) * log10junk);
junk3 = thick .* hydrconjunk;
denom = sum(junk3);
hyd_trans_profile = zeros(istep,1);
for ijunk = 1:1:itime
%    hyd_trans_profile(ijunk) = 0.5 * (1 + sum(junk3(ijunk:end)) / denom);   % hydraulic transmissivity profile
    hyd_trans_profile(ijunk) = params.Roff_c3 + (1-params.Roff_c3) * sum(junk3(ijunk:end)) / denom;   % hydraulic transmissivity profile
end

wfps_c1a = 0.03;
wfps_c2a = 0.5;
wfps_c3a = 20;
zstar = wfps_c1 * onevec +(wfps_c2a - wfps_c1a)*((dens - params.min_bulk_dens)./(dens - (params.min_bulk_dens - wfps_c3a)));
sp_yld_profile = onevec - zstar + zstar .* ((onevec - zstar) / 0.01) .* exp(-max(0.5,depth)./zstar) .*  (onevec - exp(0.01*(onevec./zstar)));
sp_yld_profile =  max(zerovec,sp_yld_profile);   % specific yield profile

M_array = M * ones(1,num_veg);
mfrac = m ./ M_array;
% 
% if (params.tf_old_new > 0.5)
%     annNPPmoss(:,1) = annNPP(:,1) + annNPP(:,6);  % modified for TOOLIK
% else
%     annNPPmoss(:,1) = annNPP(:,1);  % modified for TOOLIK
% end

annNPPmoss = sum(annNPP .* (onevec * params.mosses), 2); 
annNPPvasc = sum(annNPP,2) - annNPPmoss;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% WRITE SUMMARY RESULTS TO SCREEN
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

disp(sprintf('total age (y): %d   total mass (kg C/m2): %d   total depth (m): %d',num_years, M_TOTAL/2, Z_TOTAL));
disp(sprintf('total dC/dt (kg C/m2): %d ',M_TOTAL2));


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% WRITE OUT OUTPUT FILES: core profile, carbon time series, water time series, params, workspace
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% write out binary ('.mat') file of soil nodal temps (n-years, 12 months, 114 nodes)
%   this file can be read and plotted with 'read_plot_soiltemps_1.m'

save('soil_node_temp_month', 'soil_node_temp_month_save')

% ***************************
% conservation of mass tests
% ***************************

% j3 = annAGMASSIN + annROOTIN - annRESP*2 - del_M_tot;
% j4 = tot_npp - annAGMASSIN - annROOTIN + j5;
% results_5 = [time del_peatwater net_water_in (del_peatwater-net_water_in) j4 tot_npp annAGMASSIN annROOTIN annRESP*2 del_M_tot j3];  

% fname5 = [params.outname, '_o_conservation_test.txt'];
% fid5 = fopen(fname5,'w');  

% fprintf(fid5,'HPM6 output - conservation tests - units: water - m depth; NPP/mass - kg/m2/y \n');
% fprintf(fid5,' sim_yr del_peatwater net_water_in del-net_water tot-AG-BG  tot_NPP    AG_NPP    BG_NPP    tot_RESP  del_peat  net_of_last_4 \n');
% fprintf(fid5,'%7.1f  %10.6f  %10.6f  %11.7f %9.5f  %9.5f  %9.5f %9.5f %9.5f %9.5f %9.5f \n', results_5');
% status = fclose(fid5);

% ***************************
% final core profile
% ***************************

results_1 = [time depth M M_0 k_mean dens m del_c14_M];  

fname1 = [params.outname, '_core.txt'];
fname1h = [params.outname, '_core_header.txt'];
fid1h = fopen(fname1h,'w');  % header for profile (core) of final state

fprintf(fid1h,'HPM20 output - core of final state - units: depth & thickness: m, mass: kg/m2 or kg/m3, time: y; decomp: 1/y; WFPS: m3/m3 \n');
fprintf(fid1h,' cohort_age coh_depth coh_mass coh_m0 coh_k_mean coh_bulk_dens  m_each_PFT_for_%g_PFTs  del_c14_M \n', params.num_veg);
% fprintf(fid1,'%7.1f  %8.4f  %8.4f  %8.4f  %10.7f  %8.3f   %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f %10.3f \n', results_1');
% fprintf(fid1, results_1');
status = fclose(fid1h);
dlmwrite(fname1, results_1,'precision','%10.7f');

% results_6 = [time depth M m mfrac];  

% fname6 = [params.outname, '_o_core_by_PFT.txt'];
% fid6 = fopen(fname6,'w');  % profile (core) of final state by PFT

% fprintf(fid6,'HPM6 output - core by PFT of final state - units: depth & thickness: m, mass: kg/m2 or kg/m3, time: y; decomp: 1/y; WFPS: m3/m3 \n');
% fprintf(fid6,' cohort_age_(y) cohort_depth_(m) cohort_mass_(kg/m2) m_grass m_minhrb m_minsdg m_decshb m_brnmoss m_holsphag m_lawnsphag m_humsphag m_feath m_ombhrb m_ombsdg m_ombshb mfrac_grass mfrac_minhrb mfrac_minsdg mfrac_decshb mfrac_brnmoss mfrac_holsphag mfrac_lawnsphag mfrac_humsphag mfrac_feath mfrac_ombhrb mfrac_ombsdg mfrac_ombshb \n');
% fprintf(fid6,'%7.1f  %8.4f  %8.4f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f  %10.7f \n', results_6');
% status = fclose(fid6);

% ***************************
%  carbon time series
% ***************************

results_2 = [time annNPP annROOTIN annRESP del_c14_annRESP del_C_del_t WATER(:,7) annZ_total annM_total annWTD];

fname2 = [params.outname, '_carbon.txt'];
fname2h = [params.outname, '_carbon_header.txt'];
% fid2 = fopen(fname2,'w');  % time series of carbon dynamics
fid2h = fopen(fname2h,'w');  % header for time series of carbon dynamics

fprintf(fid2h,'HPM20 output - time series of carbon dynamics - units: depth/thickness: m; NPP: kg/m2/y; time: y \n');
fprintf(fid2h,' time     npp_each_PFT_for_%g_PFTs   root_input  del_c14_resp  ann_resp_(kgC/m2/y) ann_dC/dt_(kgC/m2/y)  delPeatHt_(cm)  peat_depth  peat_mass  WTD\n', params.num_veg);
% fprintf(fid2,'%7.1f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %10.7f %8.3f %8.3f %8.3f %8.3f \n', results_2');
%fprintf(fid2, results_2');
status = fclose(fid2h);
dlmwrite(fname2, results_2,'precision','%10.5f');

% *********** moss fraction *********************************************
% (S. Frolking)
% if (flag_bins >0)
%    fname5 = [params.outname, '_moss_fraction.txt'];
%    save(fname5, 'bin_moss_frac', '-ascii', '-tabs');
% end

% fid5 = fopen(fname5,'w');  % time series of moss fraction binned profile
% fprintf(fid2,'HPM10 output - time series of moss fraction - dimensionless \n');
% fprintf(fid2,' time     npp_grass minherb   minsedge  decidshrub brownmoss holsphag  lawnsphag humsphag  feather   ombherb   ombsedge  ombshrub root_input ann_resp_(kgC/m2/y) ann_dC/dt_(kgC/m2/y) delPeatHt_(cm)  peat_depth peat_mass WTD\n');
% fprintf(fid2,bin_moss_frac);
% status = fclose(fid5);

% ********************************************************


% ***************************
% water time series
% ***************************

results_3 = [time WATER(:,1:6) del_peat_height*100];

fname3 = [params.outname, '_water.txt'];
fid3 = fopen(fname3,'w');  % time series of water dynamics

fprintf(fid3,'HPM20 output - time series of water dynamics - units: depth/thickness: m or m/y; time: y \n');
fprintf(fid3,' time AnnDelWat(cm) annPPT   annET    annRUNOFF annRUNON   annWTD AnnDelPtHt(cm) \n');
fprintf(fid3,'%7.1f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f \n', results_3');
status = fclose(fid3);

% ***************************
% basic annual output
% ***************************
annNPP_total = sum(annNPP,2) / 2;   % divide by to to go from kg to kg C
% annNPP_herb = (annNPP(:,1) + annNPP(:,2) + annNPP(:,3) + annNPP(:,10) + annNPP(:,11)) / 2;
% annNPP_woody = (annNPP(:,4) + annNPP(:,12) + annNPP(:,13)) / 2;
% annNPP_moss = (annNPP(:,5) + annNPP(:,6) + annNPP(:,7) + annNPP(:,8) + annNPP(:,9)) / 2;

% annNPP_herb = sum(annNPP .* (onevec * (params.vasculars - params.woody)), 2) / 2;
% annNPP_woody = sum(annNPP .* (onevec * params.woody), 2) / 2;
% annNPP_moss = sum(annNPP .* (onevec * params.mosses), 2) / 2;
annNPP_herb = sum(annNPP .* (onevec * (params.vasculars - params.woody)), 2) / 2;
annNPP_woody = sum(annNPP .* (onevec * params.woody), 2) / 2;
annNPP_moss = sum(annNPP .* (onevec * params.mosses), 2) / 2;

results_13 = [time annTEMP annPPT WATER(:,3) WATER(:,4) annWTD annNPP_total annRESP del_C_del_t annNPP_herb annNPP_woody annNPP_moss annM_total/2 annZ_total];
array_dim = size(results_13);

results_13_smooth = zeros(array_dim);
results_13_smooth(:,1) = results_13(:,1);
for (jj = 2:1:array_dim(2))
%     results_13_smooth(:,jj) = smooth(results_13(:,jj),11,'loess');
    results_13_smooth(:,jj) = smooth(results_13(:,jj),4*25,'loess');

end

if (params.tf_old_new > 0.5) 
    
    moss_old = sum(sum((onevec * (params.mosses .* [old_new_ones old_new_zeros])) .* m,2));
    moss_new = sum(sum((onevec * (params.mosses .* [old_new_zeros old_new_ones])) .* m,2));
    vascular_old = sum(sum((onevec * (params.vasculars .* [old_new_ones old_new_zeros])) .* m,2));
    vascular_new = sum(sum((onevec * (params.vasculars .* [old_new_zeros old_new_ones])) .* m,2));
    
    results_13_old_new = [time annTEMP annPPT WATER(:,3) WATER(:,4) annWTD ann_NPP_old ann_NPP_new ann_resp_old ann_resp_new ann_M_old ann_M_new];
    array_dim = size(results_13_old_new);
    
    results_13_old_new_smooth = zeros(array_dim);
    results_13_old_new_smooth(:,1) = results_13_old_new(:,1);
    for (jj = 2:1:array_dim(2))
    %     results_13_old_new_smooth(:,jj) = smooth(results_13(:,jj),11,'loess');
        results_13_old_new_smooth(:,jj) = smooth(results_13_old_new(:,jj),4*25,'loess');
    
    end
end

fname13 = [params.outname, '_basic_annual_output.txt'];
fid13 = fopen(fname13,'w');  % time series of basic annual output

fprintf(fid13,'HPM20_soiltemp_output_time_series_units_DEPTH_m_MASS_kgC_m2_y_TIME_y \n');
fprintf(fid13,'  time  annTEMP   annPPT    annET   annRUNOFF  annWTD   annNPP    annDECOMP  annNCB    annNPP_herb annNPP_woody annNPP_moss peat_mass peat_height \n');
fprintf(fid13,'%7.1f %7.1f %9.5f %9.5f %9.5f %9.3f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.2f %9.5f \n', results_13');
status = fclose(fid13);

fname14 = [params.outname, '_basic_annual_smooth_output.txt'];
fid14 = fopen(fname14,'w');  % time series of basic annual output

fprintf(fid14,'HPM20_soiltemp_output_time_series_units_DEPTH_m_MASS_kgC_m2_y_TIME_y \n');
fprintf(fid14,'  time  annTEMP   annPPT    annET   annRUNOFF  annWTD   annNPP    annDECOMP  annNCB    annNPP_herb annNPP_woody annNPP_moss peat_mass peat_height \n');
fprintf(fid14,'%7.1f %7.1f %9.5f %9.5f %9.5f %9.3f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f %9.2f %9.5f \n', results_13_smooth');
status = fclose(fid14);

fname15 = [params.outname, '_basic_annual_output_old_new.txt'];
fid15 = fopen(fname15,'w');  % time series of basic annual output

if (params.tf_old_new > 0.5)
    fprintf(fid15,'HPM20_soiltemp_output_time_series_units_DEPTH_m_MASS_kgC_m2_y_TIME_y \n');
    fprintf(fid15,'  time  annTEMP   annPPT    annET   annRUNOFF  annWTD   annNPP_old annNPP_new  annRESP_old  annRESP_new  ann_M_old  ann_M_new \n');
    fprintf(fid15,'%7.1f %7.1f %9.5f %9.5f %9.5f %9.3f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f \n', results_13_old_new');
    status = fclose(fid15);

    fname16 = [params.outname, '_basic_annual_output_old_new_smooth.txt'];
    fid16 = fopen(fname16,'w');  % time series of basic annual output

    fprintf(fid16,'HPM20_soiltemp_output_time_series_units_DEPTH_m_MASS_kgC_m2_y_TIME_y \n');
    fprintf(fid16,'  time  annTEMP   annPPT    annET   annRUNOFF  annWTD   annNPP_old annNPP_new  annRESP_old  annRESP_new  ann_M_old  ann_M_new \n');
    fprintf(fid16,'%7.1f %7.1f %9.5f %9.5f %9.5f %9.3f %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f \n', results_13_old_new_smooth');
    status = fclose(fid16);

end


% ***************************
% run parameters
% ***************************

fname4 = [params.outname, '_params.txt'];
fid4 = fopen(fname4,'w');  % parameters

fprintf(fid4,'HPM20 output - parameters - units: depth/thickness: m, mass: kg/m2 or kg/m3, time: y; decomp: 1/y \n');
% fprintf(fid4,params); 
fprintf(fid4,'output file name           %8s   \n',params.outname);
fprintf(fid4,'climate input file name    %8s   \n',params.outname);
fprintf(fid4,'simulation length [y]      %g   \n',num_years);
fprintf(fid4,'simulation start [yr BP]   %g   \n',params.sim_start);
fprintf(fid4,'simulation end [yr BP]     %g   \n',params.sim_end);
fprintf(fid4,'GIPL flag [1=yes soil T]   %g   \n',params.gipl_flag);
fprintf(fid4,'GFDL flag [1=CM3, 2=ESM2M] %g   \n',params.gfdl_model_flag);
fprintf(fid4,'start sim depth [m]        %6.2f \n',params.start_depth);
fprintf(fid4,'14C decay [e-folding y]    %6.2f \n',params.tau_c14);

fprintf(fid4,'ann_temp [C]               %6.2f \n',params.ann_temp);
fprintf(fid4,'ann_ppt [m/y]              %6.2f \n',params.ann_ppt);
fprintf(fid4,'initialization WTD [m]     %6.3f \n',params.wtd_0);
fprintf(fid4,'initialization PD [m]      %6.3f \n',params.start_depth);
fprintf(fid4,'ET_0    [m/y]              %6.2f \n',params.ET_0);
fprintf(fid4,'Roff_c1                    %6.2f \n',params.Roff_c1);
fprintf(fid4,'Roff_c2                    %6.2f \n',params.Roff_c2);
fprintf(fid4,'Roff_c2a                   %6.2f \n',params.Roff_c2a);
fprintf(fid4,'Roff_c3                    %6.2f \n',params.Roff_c3);
fprintf(fid4,'Roff_c4                    %6.2f \n',params.Roff_c4);
fprintf(fid4,'runon_c1                   %6.2f \n',params.runon_c1);
fprintf(fid4,'runon_c2                   %6.2f \n',params.runon_c2);
fprintf(fid4,'runon_c3                   %6.2f \n',params.runon_c3);
fprintf(fid4,'ET_wtd_1                   %6.2f \n',params.ET_wtd_1);
fprintf(fid4,'ET_wtd_2                   %6.2f \n',params.ET_wtd_2);
fprintf(fid4,'ET_wtd_3                   %6.2f \n',params.ET_wtd_3);
fprintf(fid4,'ET_param                   %6.2f \n',params.ET_param);
fprintf(fid4,'rootin_c3                  %6.2f \n',params.rootin_c3);
fprintf(fid4,'rootin_c4                  %6.2f \n',params.rootin_c4);
fprintf(fid4,'rootin_c5                  %6.2f \n',params.rootin_c5);
fprintf(fid4,'rootin_alpha               %6.2f \n',params.rootin_alpha);
fprintf(fid4,'rootin_d80                 %6.2f \n',params.rootin_d80);
fprintf(fid4,'wfps_c1                    %6.2f \n',params.wfps_c1);
fprintf(fid4,'wfps_c2                    %6.2f \n',params.wfps_c2);
fprintf(fid4,'wfps_c3                    %6.2f \n',params.wfps_c3);
fprintf(fid4,'wfps_opt                   %6.2f \n',params.wfps_opt);
fprintf(fid4,'wfps_sat_rate              %6.2f \n',params.wfps_sat_rate);
fprintf(fid4,'wfps_min_rate              %8.4f \n',params.wfps_min_rate);
fprintf(fid4,'wfps_curve                 %6.2f \n',params.wfps_curve);
fprintf(fid4,'dens_c1                    %6.2f \n',params.dens_c1);
fprintf(fid4,'dens_c2                    %6.2f \n',params.dens_c2);
fprintf(fid4,'min_bulk_dens [kg/m3]      %6.2f \n',params.min_bulk_dens);
fprintf(fid4,'del_bulk_dens [kg/m3]      %6.2f \n',params.del_bulk_dens);
fprintf(fid4,'OM_bulk_dens  [kg/m3]      %6.2f \n',params.OM_dens);
fprintf(fid4,'anoxic_scale_length [m]    %6.2f \n \n',params.anoxia_scale_length);
fprintf(fid4,'num_veg    %6.2f \n',params.num_veg);
fprintf(fid4,'lag years for vascular WTD %6.2f \n',params.lag_years);
fprintf(fid4,'max_total_NPP[kg/m2/y]    %6.2f \n',params.max_npp);
fprintf(fid4,'NPPQ10 value              %6.3f \n',params.NPP_Q10);
fprintf(fid4,'                      grs    minh   mins   mnshr  wtms   hols   lawn   hums   fthr   ombs   ombh   ombshr   tree\n');
if (params.tf_old_new > 0.5) 
    fprintf(fid4,'NPP_relative        %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.NPP_rel(1:num_veg/2));
    fprintf(fid4,'ag_frac_npp         %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.ag_frac_npp(1:num_veg/2)');
    fprintf(fid4,'bg_frac_npp         %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.bg_frac_npp(1:num_veg/2)');
    fprintf(fid4,'WTD_opt             %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.WTD_opt(1:num_veg/2)');
    fprintf(fid4,'WTD_range_shallow   %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.WTD_range(1,(1:num_veg/2)));
    fprintf(fid4,'WTD_range deep      %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.WTD_range(2,(1:num_veg/2)));
    fprintf(fid4,'PD_opt              %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.PD_opt(1:num_veg/2));
    fprintf(fid4,'PD_opt shallow      %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.PD_range(1,(1:num_veg/2)));
    fprintf(fid4,'PD_opt deep         %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.PD_range(2,(1:num_veg/2)));
    fprintf(fid4,'decomp k_0          %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.k_0(1:num_veg/2));
else
    fprintf(fid4,'NPP_relative        %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.NPP_rel);
    fprintf(fid4,'ag_frac_npp         %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.ag_frac_npp');
    fprintf(fid4,'bg_frac_npp         %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.bg_frac_npp');
    fprintf(fid4,'WTD_opt             %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.WTD_opt');
    fprintf(fid4,'WTD_range_shallow   %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.WTD_range(1,:));
    fprintf(fid4,'WTD_range deep      %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.WTD_range(2,:));
    fprintf(fid4,'PD_opt              %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.PD_opt);
    fprintf(fid4,'PD_opt shallow      %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.PD_range(1,:));
    fprintf(fid4,'PD_opt deep         %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.PD_range(2,:));
    fprintf(fid4,'decomp k_0          %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f \n', params.k_0);   
end
fprintf(fid4,'\n total age (y): %d   total mass (kg C/m2): %d   total depth (m): %d',num_years, M_TOTAL/2, Z_TOTAL);
%fprintf(fid4,'\n total dC/dt (kg C/m2): %d ',M_TOTAL2);

status = fclose(fid4);

% ***************************
% run workspace variable arrays as '.mat' file
% ***************************

workspace_filename = [params.outname, '_ws'];
save(workspace_filename);  
    
% **************************
% some other arrays
% ***************************

fname5 = [params.outname, '_precip_used.csv'];
dlmwrite(fname5,REC_annppt);

fname6 = [params.outname, '_age_depth.csv'];
dlmwrite(fname6,age_depth);

fname62 = [params.outname, '_age_depth2.csv'];
dlmwrite(fname62,age_depth2);

fname7 = [params.outname, '_old_new_profiles.csv'];
dlmwrite(fname7,profiles);

% --- assemble final 300-year monthly layer temperature time series

soil_layer_monthly_temps_300year = zeros(300*12,63);
for jyear = 1:1:300
    for jmonth = 1:1:12
        soil_layer_monthly_temps_300year((jyear-1)*12 + jmonth, :) = 0.5 *...
            (soil_node_temp_month_save(end-300+jyear,jmonth,1:63) + ...
             soil_node_temp_month_save(end-300+jyear,jmonth,2:64));
    end
end

% fname63 = [params.outname, '_soil_layer_temps.csv'];
% dlmwrite(fname63,soil_layer_temp_month_save);

if (flag_bins > 0.5 && flag_bins < 1.5)
    fname8 = [params.outname, '_bin_M_star.csv'];
    dlmwrite(fname8,bin_M_star);
end
if (flag_bins > 1.5)
    fname7 = [params.outname, '_bin_moss_frac.csv'];
    dlmwrite(fname7,bin_moss_frac);
end

% dlmwrite('my_data.out',A, ';')
% PFT 1: grasses
% PFT 2: minerotrophic herb
% PFT 3: minerotropic sedge
% PFT 4: deciduous shrub (woody?)
% PFT 5: minerotrophic wet moss (non-sphagnum)
% PFT 6: hollow sphagnum
% PFT 7: lawn sphagnum
% PFT 8: hummock sphagnum
% PFT 9: feathermoss (& lichen?)
% PFT 10: ombrotrophic herb
% PFT 11: ombrotrophic sedge
% PFT 12: evergreen shrub (woody?)
% PFT 13: trees

