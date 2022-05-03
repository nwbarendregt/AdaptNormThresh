% Posterior_Sample_UGM.m
% Performs MCMC sampling used for fitting UGM model to tokens task data 
% from Barendregt et al., 2022.

clear

% Define chain length parameters and thinning frequency:
N_samples = 1e4; N_chain = 10; N_thin = 5;

% Define prior to enforce non-negativity of parameters:
thresh_min = 0; thresh_max = 25;
a_min = 0; a_max = 5;
sigma_min = 0; sigma_max = 5;
tau_min = 0.1; tau_max = 0.5;
mn_min = 0; mn_max = 5;
prior = [thresh_min a_min sigma_min tau_min mn_min];

% Define task parameters to simulate model:
Nt = 15; t_d = 0.170;
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
post_samp = NaN(N_samples/N_thin,5,N_chain);

% Define covariance matrix of Gaussian proposal distribution (found 
% experimentally by tuning):
S = 0.1*[10 0 0 0 0; 0 0.5 0 0 0; 0 0 0.1 0 0; 0 0 0 0.001 0; 0 0 0 0 0.05];

% Load data from burn-in:
load('model_fit_UGM.mat');

for i = 1:N_chain
    post_samp_i = NaN(N_samples/N_thin,5);

    % Seed chains with final states and likelihoods from burn-in:
    theta = model_fit_UGM(sub_ind,speed_ind).chain_init(:,i);
    L = model_fit_UGM(sub_ind,speed_ind).L_init(i);
    
    % Evolve chain using MCMC, thinning chain as defined above:
    for j = 1:(N_samples/N_thin)
        for k = 1:N_thin
            [theta,L] = MCMC_UGM(theta,L,Nt,Sub_T,stim,S,prior);
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
% load('model_fit_UGM.mat')
% model_fit_UGM(sub_ind,speed_ind).samples = post_samp;
% save('model_fit_UGM.mat','model_fit_UGM');