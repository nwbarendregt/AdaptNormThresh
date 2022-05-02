% RC_Noise_Strength_Sweep.m
% Calculates model reward rates for reward change task with varied noise
% strengths from Barendregt et al., 2022

clear

% Load noise-free reward rate data from reward change task:
load('RC_Const_Opt_Data.mat','Const_Opt_theta')
load('RC_UGM_Opt_Data.mat','UGM_Opt_theta')

% Define simulation parameters for reward change task:
T = 5; dt = 0.005;  t_i = 1;
dg = 0.001;
m = 5; c = 1;
R_1 = linspace(1,10); R_2 = 11-R_1;
R_ind = [23 78];

% Define limits on noise source strengths:
sigma_min = 0; sigma_max = 5;
mn_min = 0; mn_max = 0.25;
sigma = linspace(sigma_min,sigma_max);
mn = linspace(mn_min,mn_max);

% Define simulation parameters for empirical reward rates:
N_trial = 1e4;

% Pre-allocate reward rate storage for varied reward changes and noise
% levels:
NB_RR = NaN(length(sigma),length(mn),length(R_ind));
Const_RR = NaN(length(sigma),length(mn),length(R_ind));
UGM_RR = NaN(length(sigma),length(mn),length(R_ind));

for k = 1:length(R_ind)

    % Construct reward timeseries:
    R = NaN(1,T/dt+1); R(1:100) = R_1(R_ind(k)); R(101:end) = R_2(R_ind(k));

    % Construct optimal models for given reward timeseries:
    [NB_thresh,~] = RC_Bellmans(T,dt,t_i,dg,m,@(t) c,R);
    Const_thresh = Const_Opt_theta(R_ind(k))*ones(1,T/5/dt+1);
    UGM_theta = UGM_Opt_theta(:,R_ind(k));

    for i = 1:length(sigma)
        for j = 1:length(mn)

            % Calculate trial data for synthetic subject:
            [y,p] = RDMD_trial_generate(m,T/5,dt,sigma(i),N_trial);

            % Pre-allocate RT, reward, and time cost storage for synethic
            % subject:
            RT = NaN(N_trial,3);
            reward = NaN(N_trial,3);
            cost = NaN(N_trial,3);
            for n = 1:N_trial

                % Perform trial block using NB model:
                [RT(n,1),C] = RDMD_sim_norm(y(n,:),T/5,dt,NB_thresh,mn(j));
                reward(n,1) = C*(R_1(R_ind(k))*(RT(n,1) < 0.5)+R_2(R_ind(k))*(RT(n,1) >= 0.5));
                cost(n,1) = c*RT(n,1);
                
                % Perform trial block using Const model:
                [RT(n,2),C] = RDMD_sim_norm(y(n,:),T/5,dt,Const_thresh,mn(j));
                reward(n,2) = C*(R_1(R_ind(k))*(RT(n,2) < 0.5)+R_2(R_ind(k))*(RT(n,2) >= 0.5));
                cost(n,2) = c*RT(n,2);
                
                % Perform trial block using UGM:
                [RT(n,3),C] = RDMD_sim_UGM(p(n,:),T/5,dt,UGM_theta(1),UGM_theta(2),sigma(i),UGM_theta(4),mn(j));
                reward(n,3) = C*(R_1(R_ind(k))*(RT(n,3) < 0.5)+R_2(R_ind(k))*(RT(n,3) >= 0.5));
                cost(n,3) = c*RT(n,3);
            end

            % Calculate and store reward rates for each model:
            RR = (mean(reward,1)-mean(cost,1))./(mean(RT,1)+t_i);
            NB_RR(i,j,k) = RR(1); Const_RR(i,j,k) = RR(2); UGM_RR(i,j,k) = RR(3);
        end
    end
end
save('RC_Noise_Sweep_Data.mat','sigma','mn','NB_RR','Const_RR','UGM_RR')