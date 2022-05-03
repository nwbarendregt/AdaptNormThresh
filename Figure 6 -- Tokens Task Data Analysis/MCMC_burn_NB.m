% MCMC_burn_NB.m
% Performs MCMC burn-in used for fitting NB model to tokens task data from
% Barendregt et al., 2022.

clear

% Define chain length parameters:
N_burn = 1e4; N_chain = 10;

% Define prior for seeding chains:
R_c_min = 0; R_c_max = 10;
c_min = 0; c_max = 1;
sigma_min = 0; sigma_max = 5;
mn_min = 0; mn_max = 5;

% Define prior to enforce non-negativity of parameters:
prior = [R_c_min c_min sigma_min mn_min];

% Define task parameters to simulate model:
Nt = 15; t_d = 0.170; R_i = -1; tol = 1e-5;
speed = t_d*1000; speed_ind = 1; % 1 for slow task, 2 for fast task.

% Load subject data:
load('trials.mat'); sub_ind = 1; % Determines which subject's data to analyze.
idSubject = [6 7 9 10 11 12 13 14 15 16 17 18 20 21 22 23 24 25 26 27];
Sub_ID = idSubject(sub_ind);
Sub_T = trials.nDecisionToken((trials.nSpeedFast == speed) & (trials.idSubject == Sub_ID));
Sub_Data = histcounts(Sub_T,-0.5:1:(Nt+0.5),'normalization','probability');

% Load and format stimulus data:
Sub_stim = trials.sTokenDirs((trials.nSpeedFast == speed) & (trials.idSubject == Sub_ID));
stim = NaN(length(Sub_stim),Nt);
for i = 1:length(Sub_stim)
    stim(i,:) = str2num(strtrim(regexprep(Sub_stim{i},'.{1}','$0 ')));
    stim(i,:) = 2*(stim(i,:)-1)-1;
end

% Pre-allocate storage for final point in each chain and likelihood of
% that final point:
chain_init = NaN(4,N_chain); L_init = NaN(1,N_chain);

% Define covariance matrix of Gaussian proposal distribution (found 
% experimentally by tuning):
S(:,:,1) = 0.9*[0.75 0 0 0; 0 0.05 0 0; 0 0 0.1 0; 0 0 0 0.01];
S(:,:,2) = 0.5*[0.1 0 0 0; 0 0.05 0 0; 0 0 0.1 0; 0 0 0 0.01];
S = S(:,:,speed_ind);

for i = 1:N_chain

    % Seed each chain using uniform prior over parameters:
    theta = [rand*(R_c_max-R_c_min)+R_c_min;...
        rand*(c_max-c_min)+c_min; rand*(sigma_max-sigma_min)+sigma_min;...
        rand*(mn_max-mn_min)+mn_min];

    % Construct model from drawn parameters:
    thresh_g = tok_Bellmans_g(Nt,t_d,theta(1),R_i,@(t) theta(2),tol);

    % Pre-allocate storage of average synthetic data generated from 
    % sampled model:
    Synth_Data = zeros(1,Nt+1);
    for n = 1:50

        % Pre-allocate storage of synthetic data generated from sampled 
        % model: 
        T_Data = NaN(1,length(Sub_T));

        % Generage synthetic response data using sampled model and subject
        % stimulus:
        for l = 1:length(Sub_T)
            [T,~] = tok_sim_norm(Nt,thresh_g,theta(3),stim(l,:));
            T_Data(l) = round(theta(4)*randn+T);
            while (T_Data(l) > Nt) || (T_Data(l) < 0)
                T_Data(l) = round(theta(4)*randn+T);
            end
        end

        % Average synthetic data over many (50, found experimentally by 
        % tuning) realizations:
        Synth_Data = Synth_Data+histcounts(T_Data,-0.5:1:(Nt+0.5),'normalization','probability');
    end

    % Add small non-zero entries to compute likelihood of sampled model:
    Synth_Data(Synth_Data == 0) = eps; Synth_Data = Synth_Data/sum(Synth_Data,'all');

    % Compute likelihood of sampled model/parameters:
    L = sum(log(Synth_Data(Sub_T+1)));

    % Evolve chain using MCMC:
    for j = 1:N_burn
        [theta,L] = MCMC_NB(theta,L,Nt,t_d,R_i,tol,Sub_T,stim,S,prior);
    end

    % Store final state of chain and likelihood of the final state:
    chain_init(:,i) = theta; L_init(:,i) = L;
end

% Save data from burn-in to use for sampling (using Posterior_Sample_NB.m):
% load('model_fit_NB.mat')
% model_fit_NB(sub_ind,speed_ind).idSubject = Sub_ID;
% model_fit_NB(sub_ind,speed_ind).speed = speed;
% model_fit_NB(sub_ind,speed_ind).chain_init = chain_init;
% model_fit_NB(sub_ind,speed_ind).L_init = L_init;
% save('model_fit_NB.mat','model_fit_NB');