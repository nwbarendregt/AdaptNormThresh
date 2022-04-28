% IRC_Threshold_Dynamics.m
% Script that generates realizations of thresholds for inferred reward 
% change task from Barendregt et al., 2022. 

clear

% Define simulation parameters:
T = 10; dt = 0.005; t = 0:dt:T; t_i = 1;
dg = 0.001;
m_s = 5; m_r = [1 5 10 50]; h_r = 1;
c = @(t) 1; R_H = 10; R_L = 8;

% Define number of samples for type of threshold change (low-to-high or
% high-to-low) and window size:
N_Sample = 100; W = 0.250;

% Pre-allocate storage for windowed threshold dynamics:
thresh_HL = cell(1,length(m_r)); thresh_LH = cell(1,length(m_r));

for i = 1:length(m_r)

    % Pre-allocate storage for windowed threshold dynamics for given m_r:
    thresh_HL_i = []; thresh_LH_i = [];
    s_HL = size(thresh_HL_i); s_LH = size(thresh_LH_i);
    while (s_HL(1) < N_Sample) || (s_LH(1) < N_Sample) % Continue until desired number of samples is achieved.

        % Generate reward timeseries following two-state Markov process:
        x_r = NaN(1,length(t)); x_r(1:end) = 2*binornd(1,0.5)-1;
        for j = 2:length(t)
            if rand < dt*h_r
                x_r(j:end) = -x_r(j-1);
            end
        end

        % Calculate belief over reward state following nonlinear DDM (see
        % Methods):
        y_r = zeros(1,length(t));
        for j = 2:length(t)
            y_r(j) = y_r(j-1)+x_r(j-1)*m_r(i)*dt-2*h_r*sinh(y_r(j-1))*dt+sqrt(2*m_r(i))*randn*sqrt(dt);
        end

        % Calculate threshold using dynamic programming:
        thresh = IRC_Bellmans(y_r,h_r,T,dt,t_i,dg,m_s,c,R_H,R_L); thresh = thresh(1:(T/5/dt+1));

        % Extract windowed threshold dynamic around reward changepoints and
        % group dynamics by the type of changepoint:
        for j = (2+W/dt):(T/5/dt+1-W/dt)
            if (x_r(j-1) == 1) && (x_r(j) == -1)
                thresh_HL_i = [thresh_HL_i;thresh((j-W/dt):(j+W/dt))];
            end
            if (x_r(j-1) == -1) && (x_r(j) == 1)
                thresh_LH_i = [thresh_LH_i;thresh((j-W/dt):(j+W/dt))];
            end
        end
        s_HL = size(thresh_HL_i); s_LH = size(thresh_LH_i);
    end
    thresh_HL{i} = thresh_HL_i(1:N_Sample,:);
    thresh_LH{i} = thresh_LH_i(1:N_Sample,:);
end

% Pre-allocate storage for windowed threshold dynamics for infinite m_r:
thresh_HL_broad = []; thresh_LH_broad = [];
s_HL = size(thresh_HL_broad); s_LH = size(thresh_LH_broad);
while (s_HL(1) < N_Sample) || (s_LH(1) < N_Sample) % Continue until desired number of samples is achieved.

    % Generate reward timeseries following two-state Markov process:
    x_r = NaN(1,length(t)); x_r(1:end) = 2*binornd(1,0.5)-1;
    for j = 2:length(t)
        if rand < dt*h_r
            x_r(j:end) = -x_r(j-1);
        end
    end

    % Calculate threshold using dynamic programming:
    thresh = IRC_Bellmans_Inf(x_r,h_r,T,dt,t_i,dg,m_s,c,R_H,R_L); thresh = thresh(1:(T/5/dt+1));

    % Extract windowed threshold dynamic around reward changepoints and
    % group dynamics by the type of changepoint:
    for j = (2+W/dt):(T/5/dt+1-W/dt)
        if (x_r(j-1) == 1) && (x_r(j) == -1)
            thresh_HL_broad = [thresh_HL_broad;thresh((j-W/dt):(j+W/dt))];
        end
        if (x_r(j-1) == -1) && (x_r(j) == 1)
            thresh_LH_broad = [thresh_LH_broad;thresh((j-W/dt):(j+W/dt))];
        end
    end
    s_HL = size(thresh_HL_broad); s_LH = size(thresh_LH_broad);
    disp([s_HL(1) s_LH(1)])
end
thresh_HL_broad = thresh_HL_broad(1:N_Sample,:);
thresh_LH_broad = thresh_LH_broad(1:N_Sample,:);

% save('IRC_Threshold_Data.mat')