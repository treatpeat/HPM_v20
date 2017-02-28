function root_in = hpm_rootin20(depthvec, thickvec, params, nppvec, zwt, alt, peatheight, onevct)
% root_in = hpm_rootin12(depthvec, thickvec, params, nppvec, zwt, alt, peatheight, onevct)

% v.12 for arbitrary number of PFTs, using params.sedge and params.vascular
% to determine root inputs by PFT

% no  change from v.8

% no  change from v.5

% root profiles follow ideas of Bauer (2004) 
    %                  grs  minh mins dshr wtms hols lawn hums fthr ombh
    %                  ombs evrs trees


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

sedge_tot_root = params.bg_frac_npp .* nppvec .* params.sedges; % [0 0 nppvec(3) 0 0 0 0 0 0 0 nppvec(11) 0 0];
non_sedge_tot_root = params.bg_frac_npp .* nppvec .* (params.vasculars - params.sedges);  %[nppvec(1) nppvec(2) 0 nppvec(4) 0 0 0 0 0 nppvec(10) 0 nppvec(12) nppvec(13)];

% make all roots 'sedge'
% sedge_tot_root = params.bg_frac_npp .* [nppvec(1) nppvec(2) nppvec(3) nppvec(4) 0 0 0 0 0 nppvec(10) nppvec(11) nppvec(12)];
% non_sedge_tot_root = [0 0 0 0 0 0 0 0 0 0 0 0];

% make all roots 'non-sedge'
% sedge_tot_root = [0 0 0 0 0 0 0 0 0 0 0 0];
% non_sedge_tot_root = params.bg_frac_npp .* [nppvec(1) nppvec(2) nppvec(3) nppvec(4) 0 0 0 0 0 nppvec(10) nppvec(11) nppvec(12)];

zstar = max(zwt, params.rootin_c3); % maximum root depth for non-sedge vascular plants
if params.pf_flag > 0.5
    zstar = min(zstar, alt);
end

% SF: new routines for root input (August 2011)

%    uniform input per layer for non-sedge roots (rather than proportional to layer thickness)
%    uniform input per layer for upper range of sedge roots (depth < 'd80' from parameters (depth to 80% of root input)
%    input proportional to layer thickness below 'd80', with total of 20% from 'd80' to 2 meters

% ***SEDGE ROOTS***

input_equal_per_layer = 1;
if (input_equal_per_layer > 0.5)
    
    number_root_layers = find(depthvec > params.rootin_d80, 1,'first')-1;
%    root_frac = 1;
 
    if (isempty(number_root_layers))
         number_root_layers = find(thickvec > 0, 1,'last');
         if (isempty(number_root_layers))
             number_root_layers = 1;
         end
%         root_frac = depthvec(number_root_layers) / (zstar + eps);
    end

    upper_sedge_root_in = 0.8 / (max(1,number_root_layers));
%    non_sedge_root_in = root_frac * (onevct / (number_root_layers));
    tf_root1 = thickvec > 0;
    tf_root2 = depthvec <= params.rootin_d80;
    tf_root3 = tf_root1 .* tf_root2;
    if (isempty(tf_root3))
        tf_root3 = 0 * thickvec;
        tf_root3(1) = 1;
    end
    upper_sedge_root_in = upper_sedge_root_in .* tf_root3;
    upper_sedge_root_in = upper_sedge_root_in * min(1, peatheight/params.rootin_d80);  % adjust total to fraction of root zone that is peat

    lower_sedge_root_in = 0.2/(exp(-params.rootin_d80) - exp(-2)) * exp(-depthvec) .* thickvec;
    tf_root4 = depthvec > params.rootin_d80;
    if (isempty(tf_root4))
        tf_root4 = 0 * thickvec;
    end
    tf_root5 = depthvec <= 2;  % no sedge roots below 2 meters
    tf_root6 = tf_root4 .* tf_root5;
    lower_sedge_root_in = lower_sedge_root_in .* tf_root6;
    
    sedge_root_in = upper_sedge_root_in + lower_sedge_root_in;
    
else
    
    sedge_root_in = thickvec .* (params.rootin_alpha * exp(-params.rootin_alpha * depthvec))...
                .* (onevct - min(onevct,fix(depthvec/params.rootin_c4)));

    norm1 = sum(sedge_root_in) / min(1,(1-exp(-params.rootin_alpha * peatheight)));
            
% sedge_root_in = sedge_root_in / (sum(sedge_root_in) + eps);  % normalize total to 1.0??
    sedge_root_in = sedge_root_in / (norm1 + eps);  % normalize total to 1.0
  
end

% -- if permafrost, compress sedge root input to active layer

zbottom = cumsum(thickvec);
tf_zbot = zbottom < alt;
sedge_root_in_pf = sedge_root_in .* tf_zbot; % restrict root input to Active Layer            
if (sum(sedge_root_in_pf) > 0)
    pf_factor = 1/sum(sedge_root_in_pf);
else
    pf_factor = 1.;
end
sedge_root_in_pf = sedge_root_in_pf * pf_factor;
sedge_root_in = sedge_root_in_pf;

% ***NON-SEDGE ROOTS***

input_equal_per_layer = 1;
if (input_equal_per_layer > 0.5)
    
    number_root_layers = find(depthvec > zstar, 1,'first')-1;
%    root_frac = 1;
 
    if (isempty(number_root_layers))
         number_root_layers = find(thickvec > 0, 1,'last');
         if (isempty(number_root_layers))
             number_root_layers = 1;
         end
%         root_frac = depthvec(number_root_layers) / (zstar + eps);
    end

    non_sedge_root_in = onevct / (max(1,number_root_layers));
%    non_sedge_root_in = root_frac * (onevct / (number_root_layers));
    tf_root1 = thickvec > 0;
    tf_root2 = depthvec <= zstar;
    tf_root = tf_root1 .* tf_root2;
    if (isempty(tf_root))
        tf_root = 0 * thickvec;
        tf_root(1) = 1;
    end
    non_sedge_root_in = non_sedge_root_in .* tf_root;
    non_sedge_root_in = non_sedge_root_in * min(1, peatheight/zstar);  % adjust total to fraction of root zone that is peat

else
    
% first version (below) uses error function to get a smooth boundary,second has uniform input to zstar
% second version lost about 5% of root mass due to discretization(?), hence divided by sum...

%     non_sedge_root_in = (thickvec/zstar) .* (onevct - 0.5*(onevct + erf((depthvec - zstar*onevct)/(sqrt(2)*params.rootin_c5))));

    non_sedge_root_in = (thickvec/zstar) .* (depthvec < zstar);

    non_sedge_root_in = non_sedge_root_in / (sum((thickvec/zstar) .* (depthvec < zstar)) + eps);  % normalize total to 1.0??
    non_sedge_root_in = non_sedge_root_in * min(1, peatheight/zstar);  % adjust total to fraction of root zone that is peat
    
%     norm2 = sum((thickvec/zstar) .* (depthvec < zstar)) / min(1,(peatheight/zstar));
%     non_sedge_root_in = non_sedge_root_in / (norm2 + eps);  % normalize total to 1.0

end

root_in = sedge_root_in * sedge_tot_root + non_sedge_root_in * non_sedge_tot_root;

%j5a = [ norm1 norm2 ]

return