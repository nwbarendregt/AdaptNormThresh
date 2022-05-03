% Posterior_Sample_NB.m
% Performs MCMC sampling used for fitting NB model to tokens task data from
% Barendregt et al., 2022.

clear

% Define chain length parameters and thinning frequency:
N_samples = 1e4; N_chain = 10; N_thin = 5;

% Define prior to enforce non-negativity of parameters:
R_c_min = 0;
c_min = 0; 
sigma_min = 0; 
mn_min = 0;
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

% Pre-allocate storage for thinned chain, representing samples of model
% posterior p(theta|data):
post_samp = NaN(N_samples/N_thin,4,N_chain);

% Define covariance matrix of Gaussian proposal distribution (found 
% experimentally by tuning):
S(:,:,1) = 0.9*[0.75 0 0 0; 0 0.05 0 0; 0 0 0.1 0; 0 0 0 0.01];
S(:,:,2) = 0.5*[0.1 0 0 0; 0 0.05 0 0; 0 0 0.1 0; 0 0 0 0.01];
S = S(:,:,speed_ind);

% Load data from burn-in:
load('model_fit_NB.mat');

for i = 1:N_chain
    post_samp_i = NaN(N_samples/N_thin,4);

    % Seed chains with final states and likelihoods from burn-in:
    theta = model_fit_NB(sub_ind,speed_ind).chain_init(:,i);
    L = model_fit_NB(sub_ind,speed_ind).L_init(i);

    % Evolve chain using MCMC, thinning chain as defined above:
    for j = 1:(N_samples/N_thin)
        for k = 1:N_thin
            [theta,L] = MCMC_NB(theta,L,Nt,t_d,R_i,tol,Sub_T,stim,S,prior);
        end
        post_samp_i(j,:) = theta;
    end

    % Store thinned chain as as samples of model posterior:
    post_samp(:,:,i) = post_samp_i;
end

% Reshape thinned chains to create histogram approximation of model
% posterior:
X = [];
for i = 1:N_chain
    X = [X;post_samp(:,:,i)];
end
post_samp = X;

% Save posterior samples for secondary analysis (find MLE, model
% comparison, etc.):
% load('model_fit_NB.mat')
% model_fit_NB(sub_ind,speed_ind).samples = post_samp;
% save('model_fit_NB.mat','model_fit_NB');