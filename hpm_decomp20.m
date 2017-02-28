function wfps_fact = hpm_decomp20(depthvec, wtdnow, annwfps, params, onevct, epsvct)
% wfps_fact = hpm_decomp20(depthvec, wtdnow, annwfps, params, onevct, epsvct)
% no change from v.8

% no change from v.5

% temperature factor not yet developed (set equal to 1.0) ***SEE OTHER VERSION
% (J.Talbot)

% WFPS factor:
 
   % above WT: parabolic function, optimum at WFPS = 0.4; drops to 0.5 for WFPS = 0; drops to 0.1 for WFPS = 1.0
   %      see spreadsheet 'simple decomp modeling.xls'
   
   % below WT: exponential drop from 0.1 to 0.01 with depth (ideas of Christian Blodau)
   %   scale-depth for decline should decrease with boginess; initially set to 0.5 m
   
% inputs are all column vectors except params (input parameters)

% depthvec = vector of cohort depths [m]
% wtdnow = water table depth [m] this year.
% onevct = vector of ones
% epsvct = vector of eps (tiny positive values)

wfps_fact = (1.0 - ((annwfps - params.wfps_opt).^2)/(4 * params.wfps_curve)) .* (depthvec < wtdnow)...
            + (depthvec >= wtdnow) .* (params.wfps_min_rate + (params.wfps_sat_rate - params.wfps_min_rate)...
            * exp(-(depthvec - wtdnow) / params.anoxia_scale_length));

        % SF: modification to eliminate effect of WFPS in unsaturated zone
%        (debugging test only -- do not use)

% wfps_fact = (1.0) .* (depthvec < wtdnow)...
%            + (depthvec >= wtdnow) .* (params.wfps_min_rate + (params.wfps_sat_rate - params.wfps_min_rate)...
%            * exp(-(depthvec - wtdnow) / params.anoxia_scale_length));

%***NO TEMPERATURE EFFECT ON DECOMPOSITION YET***
%  temp_fact = onevct;

return