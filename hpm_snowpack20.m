function temperature = hpm_snowpack10(istep, params)

% function to generate an annual snowpack depth time series to drive HPM
%    initial goal is for developing a new permafrost routine

% Initially assume there is a 250-yr interval reconstruction as for
%      Mer Bleue precipitation reconstruction from Muller et al. (2003)

snowdata_0_2750BP =    [ -6. -6. -6. -6. -6. -6. -6. -6. -6. -6. -6. -6. ];
snowdata_3000_5750BP = [ -6. -6. -6. -6. -6. -6. -6. -6. -6. -6. -6. -6. ];
snowdata_6000_8500BP = [ -6. -6. -6. -6. -6. -6. -6. -6. -6. -6. -6. ];

snowconfint_0_2750BP =   [ 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 ];
snowconfint_3000_5750BP= [ 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 ];
snowconfint_6000_8500BP= [ 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 ];

snow_mean = [snowdata_0_2750BP snowdata_3000_5750BP snowdata_6000_8500BP];
snow_confint = [snowconfint_0_2750BP snowconfint_3000_5750BP snowconfint_6000_8500BP];

MB_annsnow_mean = interp1(0:250:8500, snow_mean, 1:8500, 'pchip');
MB_annsnow_ci = interp1(0:250:8500, snow_confint, 1:8500, 'pchip');

MB_annsnow_mean = fliplr(MB_annsnow_mean);
MB_annsnow_ci = fliplr(MB_annsnow_ci);

random = zeros(istep,1);
test = randn(istep,1);
% random2 = cumsum(test);
% test = 0.7 * test;
% random1 = cumsum(test) + (0.3/0.7)*test;

random(1) = test(1);
for jjj = 2:1:istep
    random(jjj) = params.ppt_rand_persist * random(jjj-1) + test(jjj);
end

m1 = max(random);
m2 = min(random);
m3 = max(abs(m1),abs(m2));
random = random / (m3+eps); % random normalized to max of 1 or min of -1 

temperature = MB_annsnow_mean + MB_annsnow_ci .* (random * params.ppt_amp1);

return

 
        