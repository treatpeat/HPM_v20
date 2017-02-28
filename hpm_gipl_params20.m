% parameters used by gipl soil temperature routine

% v.20: making version for version control - June 2015

soilLayerDepth1 = (0.0  : 0.05 : 1.0);      soilLayerDepth2 = (1.1  : 0.10 : 3.0);
soilLayerDepth3 = (3.2  : 0.20 : 6.0);      soilLayerDepth4 = (6.5  : 0.50 : 20.0);
soilLayerDepth5 = (21.0 : 1.0  : 30.0);     soilLayerDepth6 = (32.0 : 2.0  : 50.0);
soilLayerDepth7 = (55.0 : 5.0  : 100.0);

soilLayerDepth = [soilLayerDepth1 soilLayerDepth2 soilLayerDepth3 soilLayerDepth4 ...
                    soilLayerDepth5 soilLayerDepth6 soilLayerDepth7];
NumberOfSoilComputationNodes = length(soilLayerDepth);
ndepth = find(soilLayerDepth > 10., 1) -1;  % # layers above 10 m (above bedrock)

% ----------- phase change temp range: (Tfr - FIT) to (Tfr + FIT)
FIT = 0.08; % half-width of phase change temperature range (deg C)
Tfr= -FIT;  % 0.; % middle of phase change temperature range (deg C)

% ----------- Thermal conductivity (W/m/K)
Lpeat = 0.06;
Lmin = 2.0;
Lwat = 0.57;
Lice = 2.2;
Lair = 0.025;
% ----------- volumetric heat capacity (J/m3/K)
Cpeat = 0.58e6;
Cmin = 0.9e6;  % dry mineral soil (from Wisser et al.)
Cwat = 4.18e6;
Cice = 1.9e6;
Cair = 1.25e3;
% ----------- water latent heat of fusion (J/m3)
latentheat = 334e6;

apparent_heat_cap = latentheat / (2*FIT); % (J/m3/K)

Wunf = 0.02; % uses 2% for unfrozen water (Treat code)

%               jan feb mar apr may jun jul aug sep oct nov dec
snowDensity = [ 350 400 450 500 500 500 500 500 200 200 250 300 ]; % kg/m3
%  Wisser et al. 2011 starts at 150 and increases by 3 kg/m3/d  (seems too high for AK by late spring

% ------------ initial temperature profile 'T_init' at depths 'D_init' (January)
D_init = [0	0.5	1    2	3	  4	5	  10	15	  20	40	  60	80	100]; %soilLayerDepth (temperature nodes, m);
T_init = [-10 -5 -2 -1.16 -1.78 -2.18 -2.40 -2.36 -2.22 -2.09 -1.51 -0.98 -0.64 -0.64];  % TOOLIK node temperatures
% T_init = [0 3 5 6 6 6 6 6 6 6 6 6 6 6];  % MER BLEUE node temperatures
Dn_init = length(D_init);

% interpolate to GIPL2 grid (soilLayerDepth)
T_init_gipl = interp1(D_init,T_init,soilLayerDepth,'pchip');

ALFA0 = 20.14;  % snow heat transfer parameter (W/m2/K)
    
max_ABS = 1.41e-6;    % numerical iteration threshold
iter0 = 21;  % numerical iteration threshold
E0 = 1.40e-6; % numerical iteration threshold
G0 = 0.0005; % original value is 0.015; % geothermal heat flux? (W/m2)

save('hpm20_gipl_param_vals','Lpeat','Lmin','Lwat','Lice','Lair',...
    'Cpeat','Cmin','Cwat','Cice','Cair',...
    'soilLayerDepth','NumberOfSoilComputationNodes','ndepth',...
    'Tfr','FIT','latentheat','apparent_heat_cap','Wunf',...
    'snowDensity','ALFA0','max_ABS','iter0','E0','G0','T_init_gipl');
