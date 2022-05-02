% RC_Const_Opt.m
% Calculates optimal parameterization of the noise-free Const model, and 
% associated reward rate, for the reward change task from 
% Barendregt et al., 2022.

clear

% Define simulation parameters for the reward change task:
T = 5; dt = 0.005;  t_i = 1;
dg = 0.001;
m = 5; c = 1;
R_1 = linspace(1,10); R_2 = 11-R_1;

% Define simulation parameters to empirically calculate model performance:
N_trial = 1e4;

% Define model parameterization mesh to search over:
Const_thresh = linspace(0.1,5,10);

% Pre-allocate reward rate and model parameterization storage:
Const_RR_Opt = NaN(1,length(R_1));
Const_Opt_theta = NaN(1,length(R_1));

for i = 1:length(R_1)

    % Pre-allocate RT, reward, and time cost storage for model
    % parameterizations:
    RT = NaN(length(Const_thresh),N_trial);
    reward = NaN(length(Const_thresh),N_trial);
    cost = NaN(length(Const_thresh),N_trial);
    for j = 1:length(Const_thresh)

        % Construct Const model given parameterization:
        thresh = Const_thresh(j)*ones(1,T/5/dt+1);

        % Generate block data:
        y = RDMD_trial_generate(m,T/5,dt,0,N_trial);
        for n = 1:N_trial

            % Perform trial block using given Const model parameterization:
            [RT(j,n),C] = RDMD_sim_norm(y(n,:),T/5,dt,thresh,0);
            reward(j,n) = C*(R_1(i)*(RT(j,n) < 0.5)+R_2(i)*(RT(j,n) >= 0.5));
            cost(j,n) = c*RT(j,n);
        end
    end

    % Calculate maximal reward rate for given reward timeseries and store
    % optimal parameterization:
    [Const_RR_Opt(i),I] = max((mean(reward,2)-mean(cost,2))./(mean(RT,2)+t_i));
    Const_Opt_theta(i) = Const_thresh(I);
end
save('RC_Const_Opt_Data.mat','Const_RR_Opt','Const_Opt_theta','R_1','R_2');