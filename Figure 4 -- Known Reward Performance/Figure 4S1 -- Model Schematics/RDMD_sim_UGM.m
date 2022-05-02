% RDMD_sim_UGM.m
% Generates response times for the UGM with added motor noise from 
% Barendregt et al., 2022.

function [RT,C] = RDMD_sim_UGM(p,T,dt,thresh,a,sigma,tau,mn)

% Pre-allocate noise-free response time:
RT0 = NaN;

% Instantly respond with random choice if threshold is trivally-collapsed:
if thresh == 0
    RT0 = 0; C = binornd(1,0.5);
else
    t = 0; E = 0; x = 0; k = 0;

    % Determine threshold crossing time and associated choice:
    while isnan(RT0)
        if abs(x) >= thresh
            RT0 = t; C = (x > 0);
        elseif t >= T
            RT0 = T; C = (x > 0);
        else
            t = t+dt; k = k+1;
            E = E+((-E+p(k)-0.5)*dt+sigma*randn*sqrt(dt))/tau;
            x = a*E*t;
        end
    end
end
if mn>0

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