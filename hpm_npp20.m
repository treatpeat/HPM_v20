function productivity = hpm_npp20(annWTD,lagWTD,annALD,thickvec,params)
% productivity = hpm_npp20(annWTD,lagWTD,annALD,thickvec,params)

% no change from v.9
% no change from v.6

% calculates productivity of each veg type, based on WTD and parameters

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
% PFT 13: trees (J. Talbot)

% Changes in version 6 (May 2008):

% • NPP productivity response functions to WTD and PD are now skewed (in
%      some cases) -- see 'vegetation_NPP_4.m'
% • params for WTD_range and PD_range now have values for both shallower
%      and deeper sides (to generate skewness)
% • function chooses which to use based on wtdvec < or > WTD_opt, pdvec < or > PD_opt


peatdepth = sum(thickvec);
nppmin = 0.000001;       % min NPP = 1 mg/m2/y for each veg type to prevent divide by zero errors and to provide 'seed stock'

one_vec = ones(1,params.num_veg);

%             grs  minh mins dshr wtms hols lawn hums fthr ombh ombs evrs
%             trees
% WTD used for NPP is 'current' for non-vascular plants, mean of last 10 years for vascular plants

wtdvec = annWTD * params.mosses + lagWTD * params.vasculars; %+ lagWTDb * [0 0 0 0 0 0 0 0 0 0 0 0 1]; %JULIE T for wtdlagb 
pdvec  = peatdepth * one_vec;
nppminvec = nppmin * one_vec;
aldvec = annALD * one_vec;

wtd_shallow = wtdvec < params.WTD_opt;
wtd_deep = wtd_shallow == 0;

pd_shallow = pdvec < params.PD_opt;
pd_deep = pd_shallow == 0;

ald_shallow = aldvec < params.ALD_opt;
ald_deep = ald_shallow == 0;

wtd_var = wtd_shallow .* params.WTD_range(1,:) + wtd_deep .* params.WTD_range(2,:);
pd_var = pd_shallow .* params.PD_range(1,:) + pd_deep .* params.PD_range(2,:);
ald_var = ald_shallow .* params.ALD_range(1,:) + ald_deep .* params.ALD_range(2,:);

if params.pf_flag > 0.5
    NPP = params.NPP_rel .* exp(-(((wtdvec - params.WTD_opt) ./ wtd_var).^2 + ((aldvec - params.ALD_opt) ./ ald_var).^2));
else
    NPP = params.NPP_rel .* exp(-(((wtdvec - params.WTD_opt) ./ wtd_var).^2 + ((pdvec - params.PD_opt) ./ pd_var).^2));
end


productivity = max(nppminvec, NPP);        %  'kg/m2/y'

return


ald_shallow = aldvec < params.ALD_opt;
ald_deep = ald_shallow == 0;

wtd_var = wtd_shallow .* params.WTD_range(1,:) + wtd_deep .* params.WTD_range(2,:);
pd_var = pd_shallow .* params.PD_range(1,:) + pd_deep .* params.PD_range(2,:);
ald_var = ald_shallow .* params.ALD_range(1,:) + ald_deep .* params.ALD_range(2,:);

% NPP = params.NPP_rel .* exp(-(((wtdvec - params.WTD_opt) ./ wtd_var).^2 + ((pdvec - params.PD_opt) ./ pd_var).^2));
NPP = params.NPP_rel .* exp(-(((wtdvec - params.WTD_opt) ./ wtd_var).^2 + ((aldvec - params.ALD_opt) ./ ald_var).^2));

