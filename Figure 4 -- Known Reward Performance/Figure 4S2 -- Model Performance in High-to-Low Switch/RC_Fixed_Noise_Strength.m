% RC_Fixed_Noise_Strength.m

load('DR_RR_Const_Opt_data.mat','Const_Opt_theta')
load('DR_RR_UGM_Opt_data.mat','UGM_Opt_theta')

% Define simulation variables for reward change task simulations:
T = 5; dt = 0.005;  t_i = 1;
dg = 0.001;
m = 5; c = 1;
R_1 = 8; R_2 = 11-R_1; R_ind = 78;
R = NaN(1,T/dt+1); R(1:100) = R_1; R(101:end) = R_2;

% Define noise strength discretization and maximal levels of each 
% noise source:
noise_strength = 0:0.1:1;
sigma_max = 5; mn_max = 0.25;

% Define simulation parameters to construct synthetic subjects:
N_sub = 100; N_trial = 1e4;

% Pre-allocate model reward rate storage:
NB_RR = NaN(length(noise_strength),N_sub);
Const_RR = NaN(length(noise_strength),N_sub);
UGM_RR = NaN(length(noise_strength),N_sub);

% Construct optimal models given reward timeseries
[NB_thresh,~] = RC_Bellmans(T,dt,t_i,dg,m,@(t) c,R);
Const_thresh = Const_Opt_theta(R_ind)*ones(1,T/5/dt+1);
UGM_theta = UGM_Opt_theta(:,R_ind);

for i = 1:length(noise_strength)
    NB_RR_i = NaN(1,N_sub);
    Const_RR_i = NaN(1,N_sub);
    UGM_RR_i = NaN(1,N_sub);
    
    % Sample strength of each noise source given fixed total noise level:
    z = noise_strength(i)*(sigma_max+mn_max);
    sigma = (min([sigma_max,z])-max([z-mn_max,0]))*rand(1,N_sub)+max([z-mn_max,0]);
    mn = z-sigma;

    for j = 1:N_sub

        % Generate trial data for synthetic subject:
        [y,p] = RDMD_trial_generate(m,T/5,dt,sigma(j),N_trial);

        % Pre-allocate RT, reward, and time cost storage for synthetic
        % subject:
        RT = NaN(N_trial,3);
        reward = NaN(N_trial,3);
        cost = NaN(N_trial,3);
        for n = 1:N_trial

            % Perform trial block using NB model:
            [RT(n,1),C] = RDMD_sim_norm(y(n,:),T/5,dt,NB_thresh,mn(j));
            reward(n,1) = C*(R_1*(RT(n,1) < 0.5)+R_2*(RT(n,1) >= 0.5));
            cost(n,1) = c*RT(n,1);
            
            % Perform trial block using Const model:
            [RT(n,2),C] = RDMD_sim_norm(y(n,:),T/5,dt,Const_thresh,mn(j));
            reward(n,2) = C*(R_1*(RT(n,2) < 0.5)+R_2*(RT(n,2) >= 0.5));
            cost(n,2) = c*RT(n,2);
            
            % Perform trial block using UGM:
            [RT(n,3),C] = RDMD_sim_UGM(p(n,:),T/5,dt,UGM_theta(1),UGM_theta(2),UGM_theta(3)+sigma(j),UGM_theta(4),mn(j));
            reward(n,3) = C*(R_1*(RT(n,3) < 0.5)+R_2*(RT(n,3) >= 0.5));
            cost(n,3) = c*RT(n,3);
        end

        % Calculate synthetic reward rates for each model:
        RR = (mean(reward,1)-mean(cost,1))./(mean(RT,1)+t_i);
        NB_RR_i(j) = RR(1); Const_RR_i(j) = RR(2); UGM_RR_i(j) = RR(3);
    end
    NB_RR(i,:) = NB_RR_i; Const_RR(i,:) = Const_RR_i; UGM_RR(i,:) = UGM_RR_i;
end

filename = ['RC_Fixed_Noise_Strength_R1_' num2str(R_1) '.mat'];
save(filename)