% tok_sim_UGM.m
% Generates response times for tokens task using the discrete-time UGM from
% Barendregt et al., 2022.

function T = tok_sim_UGM(Nt,theta,a,sigma,tau,varargin)

% Use specified token movement stimulus if specified. Otherwise, generate
% random token movement stimulus:
if isempty(varargin)
    xi = [NaN 2*binornd(1,0.5,1,Nt)-1];
else
    xi = [NaN varargin{1}];
end

% Pre-allocate noise-free response time:
T = NaN;

% Define timestep for discretization:
n = 0; t = 0; dt = 0.200;

% Instantly respond with random choice if threshold is trivally-collapsed:
if theta == 0
    T = 0;
else
    E = 0; x = 0;

    % Determine threshold crossing time and associated choice:
    while isnan(T)
        n = n+1; t = t+dt;
        if n>Nt
            T = Nt;
        else
            p = prob_choice(sum(xi(1:n)==1),sum(xi(1:n)==-1),Nt);
            E = E+((-E+p-0.5)*dt+sigma*randn*sqrt(dt))/tau;
            x = a*E*t;
            if x>=theta
                T = n;
            elseif x<=-theta
                T = n;
            end
        end
    end
end