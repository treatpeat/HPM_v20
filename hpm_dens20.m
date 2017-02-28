function density = hpm_dens20(mass_star,mass_overlying,params,onevct)
% density = hpm_dens20(mass_star,mass_overlying,params,onevct)

% no change from v.8

% no change from v.5

% function calculates peat density at 'depth' in profile 
% uses error function to get shape

%  mass_star = fraction of original slow pool litter mass remaining
%  mass_overlying = total mass of all overlying cohorts
%  density = peat bulk density [kg/m3]
%  min_bulk_dens = surface (assumed minimum) bulk density [kg/m3]
%  del_bulk_dens = increase in bd [kg/m3] from surface to base (assumed maximum)
%  c1 = controls humification (m*) at which bulk density transition occurs [--] 
%  c2 = controls steepness of bulk density transition [--]

% NEED TO FIGURE OUT HOW TO USE 
%   1. MASS_OVERLYING???
%   2. Shrub root biomass (rigid structure--read malmer paper)???
%   3. peat water content (incompressibility)????

density = params.min_bulk_dens * onevct + params.del_bulk_dens .* (onevct - 0.5*(onevct + erf((mass_star - params.dens_c1*onevct)/params.dens_c2/sqrt(2))));

return