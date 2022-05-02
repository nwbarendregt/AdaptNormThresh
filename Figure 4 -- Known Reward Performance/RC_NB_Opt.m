% RC_NB_Opt.m
% Calculates reward rate of the noise-free NB model for the reward change
% task from Barendregt et al., 2022.

clear

% Define simulation parameters for the reward change task:
T = 5; dt = 0.005;  t_i = 1;
dg = 0.001;
m = 5; c = 1;
R_1 = linspace(1,10); R_2 = 11-R_1;

% Pre-allocate reward rate storage:
NB_RR_Opt = NaN(1,length(R_1));

for i = 1:length(R_1)

    % Construct reward timeseries:
    R = NaN(1,T/dt+1); R(1:100) = R_1(i); R(101:end) = R_2(i);

    % Calculate model reward rate using dynamic programming:
    [~,NB_RR_Opt(i)] = dynamic_reward_det(T,dt,t_i,dg,m,@(t) c,R);
end
save('RC_NB_Opt_Data.mat','NB_RR_Opt','R_1','R_2');