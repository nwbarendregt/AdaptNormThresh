% SC_Threshold_Motifs.m
% Script used to generate data for Fig. 3C in Barendregt et al., 2022.

clear

% Define simulation parameters:
T = 5; dt = 0.005;  t_i = 1;
dg = 0.001;
R = 5; c = @(t) 1;

% Define pre- and post-change rewards:
m_1 = linspace(1,10); m_2 = linspace(1,10);

% Pre-allocate threshold classification matrix and threshold storage:
type = NaN(length(m_1),length(m_2));
thresh = NaN(length(m_1),length(m_2),201);

for i = 1:length(m_1)
    type_i = NaN(1,length(m_2));
    thresh_i = NaN(length(m_2),201);
    for j = 1:length(m_2)

        % Construct SNR timeseries:
        m = NaN(1,T/dt+1); m(1:100) = m_1(i); m(101:end) = m_2(j);

        % Calculate normative thresholds using Bellman's equation:
        theta = SC_Bellmans(T,dt,t_i,dg,m,c,R); theta = theta(1:201);
        thresh_i(j,:) = theta;

        % Classify threshold motif:
        if sum(diff(theta)==0)==(length(theta)-1)
            type_i(j) = 1; % Constant Threshold
        elseif (sum(diff(theta)<=0)==(length(theta)-1))
            type_i(j) = 2; % Decreasing
        else
            type_i(j) = 3; % Increasing
        end
    end
    type(i,:) = type_i;
    thresh(i,:,:) = thresh_i;
end
save('SC_Thresholds_Motifs.mat');