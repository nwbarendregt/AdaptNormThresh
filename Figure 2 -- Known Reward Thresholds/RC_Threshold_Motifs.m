% RC_Threshold_Motifs.m
% Script used to generate data for Fig. 1C in Barendregt et al., 2022.

clear

% Define simulation parameters:
T = 5; dt = 0.005;  t_i = 1;
dg = 0.001;
m = 5; c = @(t) 1;

% Define pre- and post-change rewards:
R_1 = linspace(1,10); R_2 = linspace(1,10);

% Pre-allocate threshold classification matrix and threshold storage:
type = NaN(length(R_1),length(R_2));
thresh = NaN(length(R_1),length(R_2),T/5/dt+1);

for i = 1:length(R_1)
    type_i = NaN(1,length(R_2));
    thresh_i = NaN(length(R_2),201);
    for j = 1:length(R_2)
        
        % Construct reward timeseries:
        R = NaN(1,T/dt+1); R(1:100) = R_1(i); R(101:end) = R_2(j);

        % Calculate normative thresholds using Bellman's equation:
        theta = RC_Bellmans(T,dt,t_i,dg,m,c,R); theta = theta(1:(T/5/dt+1));
        thresh_i(j,:) = theta;

        % Classify threshold motif:
        if theta(1)==0
            type_i(j) = 1; % Motif vi (trival thresholds)
        elseif sum(diff(theta)==0)==(length(theta)-1)
            type_i(j) = 5; % Motif iii
        elseif (sum(diff(theta)<=0)==(length(theta)-1))
            type_i(j) = 2; % Motif v
        elseif (sum(theta(1:100)==Inf)==100 && sum(diff(theta(101:end))==0)==(length(theta(101:end))-1))
            type_i(j) = 3; % Motif i
        else
            type_i(j) = 4; % Motif ii and iv
        end
    end
    type(i,:) = type_i;
    thresh(i,:,:) = thresh_i;
end

save('RC_Threshold_Motifs.mat');