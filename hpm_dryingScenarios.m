% Drying scenarios, etc. 
% moved from main code, 2017-02-21

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  ADD A 'DRAINAGE DITCH' FOR FINAL 150 YEARS OF SIMULATION; 
%  JT suggests letting it decay over 250 years if not maintained
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%   if (itime > num_years - 150.5)
%       params.Roff_c2 = Roff_c2_orig + 0.15;
%   else
%       params.Roff_c2 = Roff_c2_orig;
%   end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  DRYING BY INCREASING ET_0 
%    increase by 22% over 250 years then persist for 1750 years to end of simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%    if (itime > num_years - 2000.5 && itime < num_years - 1749.5)
%        params.ET_0 = ET_0_orig * (1 + 0.33 * (itime - (num_years - 2000)) / 250);
%    end

%    if (itime > 5999.5 && itime < 6500.5)
%        params.ET_0 = ET_0_orig * (1 + 0.33 * (itime - 6000) / 500);
%    end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  DRYING BY DECREASING PPT
%    decrease by 20% over 250 years then persist for 1750 years to end of simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%    annppt_modifier = 0.;
%   if (itime < num_years - 2000.5)
%        annppt_modifier = 0.;
%    end
    
%    if (itime > num_years - 2000.5 && itime < num_years - 1749.5)
%        annppt_modifier = -0.2 * mean_REC_annppt * ((itime - (num_years - 2000)) / 250);
%    end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  DRYING BY DECREASING PPT AND INCREASING ET0
%    increase/decrease by 10% over 250 years then persist for 1750 years to end of simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%    annppt_modifier = 0.;
%   if (itime < num_years - 2000.5)
%        annppt_modifier = 0.;
%    end
    
%    if (itime > num_years - 2000.5 && itime < num_years - 1749.5)
%        annppt_modifier = -0.1 * mean_REC_annppt * ((itime - (num_years - 2000)) / 250);
%        params.ET_0 = ET_0_orig * (1 + 0.1 * (itime - (num_years - 2000)) / 250);
%    end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  INTERMITTENT DRYING BY DECREASING PPT AND INCREASING ET0
%    increase/decrease by 10% over 250 years then persist for 1750 years to end of simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%    if (itime < num_years - 4500.5)
%         annppt_modifier = 0.;
%    end
%     
%    if (itime > 4000 && itime < 4249.5)
%        annppt_modifier = -0.1 * mean_REC_annppt * ((itime - 4000) / 250);
%        params.ET_0 = ET_0_orig * (1 + 0.1 * (itime - 4000) / 250);
%    elseif (itime > 4749.5 && itime < 5000.5)
%        annppt_modifier = -0.1 * mean_REC_annppt * (1 - (itime - 4750) / 250);
%        params.ET_0 = ET_0_orig * (1 + 0.1 * (1 - (itime - 4750) / 250));
%         
%    elseif (itime > 5500.5 && itime < 5749.5)
%        annppt_modifier = -0.1 * mean_REC_annppt * ((itime - 5500) / 250);
%        params.ET_0 = ET_0_orig * (1 + 0.1 * (itime - 5500) / 250);
%    elseif (itime > 6249.5 && itime < 6500.5)
%        annppt_modifier = -0.1 * mean_REC_annppt * (1 - (itime - 6250) / 250);
%        params.ET_0 = ET_0_orig * (1 + 0.1 * (1 - (itime - 6250) / 250));
%        
%    elseif (itime > 7000.5 && itime < 7249.5)
%        annppt_modifier = -0.1 * mean_REC_annppt * ((itime - 7000) / 250);
%        params.ET_0 = ET_0_orig * (1 + 0.1 * (itime - 7000) / 250);
%    elseif (itime > 7749.5 && itime < 8000.5)
%        annppt_modifier = -0.1 * mean_REC_annppt * (1 - (itime - 7750) / 250);
%        params.ET_0 = ET_0_orig * (1 + 0.1 * (1 - (itime - 7750) / 250));
%    end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  DRYING BY INCREASING RUNOFF
%    increase over 250 years then persist for 1750 years to end of simulation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%    if (itime > num_years - 2000.5 && itime < num_years - 1749.5)
%        params.Roff_c2 = Roff_c2_orig + 0.1 * ((itime - (num_years - 2000)) / 250);
%    end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  SLOW DRYING BY INCREASING ET_0 & NO MOSS
%    increase by 25% over 2500 years (6000 to 8500)
%    moss potential NPP transfered to vascular
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%    if (itime > 5999.5)
%        params.ET_0 = ET_0_orig * (1 + 0.25 * (itime - 6000) / 2500);
%        params.NPP_rel = NPP_rel_vasc_orig * (1. + NPP_rel_orig_moss_frac);
%    end
    
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  SLOW WARMING THEN PLATEAU
%    increase by 4°C over 100 years (8200 to 8300) then hold
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%     if (itime > 8200.5)
%         if (itime < 8300.5)
%             T_mean = T_mean + 4.0 * (itime - 8200)/100.;
%         else
%             T_mean = T_mean + 4.0;
%         end
%     end
