% RDMD_sim_norm.m
% Generates response times for the NB or Const models with added motor
% noise from Barendregt et al., 2022.

function [RT,C] = RDMD_sim_norm(y,T,dt,thresh,mn)

% Pre-allocate noise-free response time:
RT0 = NaN;

% Instantly respond with random choice if threshold is trivally-collapsed:
if thresh == 0
    RT0 = 0; C = binornd(1,0.5);
else
    t = 0; k = 1;
    
    % Determine threshold crossing time and associated choice:
    while isnan(RT0)
        if t >= T
            RT0 = T; C = (y(k-1) > 0) + (y(k-1) == 0)*binornd(1,0.5);
        elseif abs(y(k)) >= thresh(k)
            RT0 = t; C = (y(k) > 0) + (y(k) == 0)*binornd(1,0.5);
        else
            k = k+1; t = t+dt;
        end
    end
end
if mn > 0

    % Filter responce time through Gaussian motor noise filter:
    RT = RT0+mn*randn;

    % Truncate Gaussian filter to ensure response in simulation domain:
    while (RT < 0) || (RT > T)
        RT = RT0+mn*randn;
    end
else
    RT = RT0;
end
end