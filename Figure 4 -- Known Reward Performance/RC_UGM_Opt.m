% RC_UGM_Opt.m
% Calculates optimal parameterization of the UGM model, and associated 
% reward rate, for the reward change task from Barendregt et al., 2022.

clear

% Define simulation parameters for the reward change task:
T = 5; dt = 0.005;  t_i = 1;
m = 5; c = 1;
R_1 = linspace(1,10); R_2 = 11-R_1;

% Define simulation parameters to empirically calculate model performance:
N_trial = 1e4;

% Define model parameterization mesh to search over:
UGM_thresh = linspace(0.1,5,10);
UGM_a = linspace(0.1,10,10);
UGM_sigma = linspace(0,5,10);
UGM_tau = logspace(-5,0,10);

% Pre-allocate reward rate and model parameterization storage:
UGM_RR_Opt = NaN(1,length(R_1));
UGM_Opt_theta = NaN(4,length(R_1));

for i = 1:length(R_1)

    % Pre-allocate reward rate storage for model parameterizations:
    RR_i = NaN(length(UGM_thresh),length(UGM_a),length(UGM_sigma),length(UGM_tau));
    for j1 = 1:length(UGM_thresh)
        for j2 = 1:length(UGM_a)
            for j3 = 1:length(UGM_sigma)
                for j4 = 1:length(UGM_tau)

                    % Pre-allocate RT, reward, and time cost storage for model
                    % parameterization:
                    RT = NaN(1,N_trial); reward = NaN(1,N_trial); cost = NaN(1,N_trial);

                    % Generate block data:
                    y = RDMD_trial_generate(m,T/5,dt,UGM_sigma(j3),N_trial);
                    for n = 1:N_trial

                        % Perform trial block using given UGM model parameterization:
                        [RT(n),C] = RDMD_sim_UGM(y(n,:),T/5,dt,UGM_thresh(j1),UGM_a(j2),UGM_tau(j4),0);
                        reward(n) = C*(R_1(i)*(RT(n) < 0.5)+R_2(i)*(RT(n) >= 0.5));
                        cost(n) = c*RT(n);
                    end
                    RR_i(j1,j2,j3,j4) = (mean(reward)-mean(cost))/(mean(RT)+t_i);
                end
            end
        end
    end

    % Calculate maximal reward rate for given reward timeseries and store
    % optimal parameterization:
    [UGM_RR_Opt(i),I] = max(RR_i,[],'all','linear');
    [th_I,a_I,s_I,t_I] = ind2sub(size(RR_i),I);
    UGM_Opt_theta(:,i) = [UGM_thresh(th_I);UGM_a(a_I);UGM_sigma(s_I);UGM_tau(t_I)];
end
save('RC_UGM_Opt_Data.mat','UGM_RR_Opt','UGM_Opt_theta','R_1','R_2');