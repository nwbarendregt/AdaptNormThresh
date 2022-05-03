% Model_Fit_NB.m
% Takes chain data from MCMC routine (Posterior_Sample_NB.m), generates the
% model posterior, and finds the maximum-likelihood model fit from
% Barendregt et al., 2022.

clear

% Load MCMC and subject data:
load('model_fit_NB.mat'); subs = 20; speeds = 2;
load('trials.mat')

% Define task parameters to simulate model:
Nt = 15; R_i = -1; tol = 1e-5;

% Define binning for empirical model posterior:
nbin = [32 32 32 32];

for speed_ind = 1:speeds
    for sub_ind = 1:subs

        % Pull task data from MCMC data:
        speed = model_fit_NB(sub_ind,speed_ind).speed; t_d = speed/1000;
        post_samp = model_fit_NB(sub_ind,speed_ind).samples;
        
        % Load subject data:
        Sub_ID = model_fit_NB(sub_ind,speed_ind).idSubject;
        Sub_T = trials.nDecisionToken((trials.nSpeedFast == speed) & (trials.idSubject == Sub_ID));
        Sub_Data = histcounts(Sub_T,-0.5:1:(Nt+0.5),'normalization','probability');
        
        % Load and format stimulus data:
        Sub_stim = trials.sTokenDirs((trials.nSpeedFast == speed) & (trials.idSubject == Sub_ID));
        stim = NaN(length(Sub_stim),Nt);
        for i = 1:length(Sub_stim)
            stim(i,:) = str2num(strtrim(regexprep(Sub_stim{i},'.{1}','$0 ')));
            stim(i,:) = 2*(stim(i,:)-1)-1;
        end
        
        % Construct empirical model posterior as multi-dimensional
        % histogram:
        edges = {linspace(0.9*min(post_samp(:,1)),1.1*max(post_samp(:,1)),nbin(1)+1),...
            linspace(0.9*min(post_samp(:,2)),1.1*max(post_samp(:,2)),nbin(2)+1),...
            linspace(0.9*min(post_samp(:,3)),1.1*max(post_samp(:,3)),nbin(3)+1),...
            linspace(0.9*min(post_samp(:,4)),1.1*max(post_samp(:,4)),nbin(4)+1)};
        [post,~,mid,~] = histcn(post_samp,edges{1},edges{2},edges{3},edges{4});
        
        % Find maximum-likelihood parameters:
        [~,I] = max(post,[],'all','linear');
        [R_c_I,c_I,sigma_I,mn_I] = ind2sub(size(post),I);
        MLE = [mid{1}(R_c_I) mid{2}(c_I) mid{3}(sigma_I) mid{4}(mn_I)];
        model_fit_NB(sub_ind,speed_ind).MLE = MLE;
        
        % Pre-allocate storage of empirical RT distribution generated from 
        % maximum-likelihood parameters:
        Fit_Data = zeros(1,Nt+1);

        % Construct maximum-likelihood model:
        thresh_g = tok_Bellmans_g(Nt,t_d,MLE(1),R_i,@(t) MLE(2),tol);
        for i = 1:50
            Fit_Data_T = NaN(1,length(Sub_T));

            % Generage synthetic response data using sampled model and subject
            % stimulus:
            for j = 1:length(Sub_T)
                [T,~] = tok_sim_norm(Nt,thresh_g,MLE(3),stim(j,:));
                Fit_Data_T(j) = round(MLE(4)*randn+T);
                while (Fit_Data_T(j) > Nt) || (Fit_Data_T(j) < 0)
                    Fit_Data_T(j) = round(MLE(4)*randn+T);
                end
            end

            % Average RT distribution over many (50, found experimentally 
            % by tuning) realizations:
            Fit_Data = Fit_Data+histcounts(Fit_Data_T,-0.5:1:(Nt+0.5));
        end

        % Add small non-zero entries to compute likelihood of sampled model:
        Fit_Data(Fit_Data == 0) = eps; Fit_Data = Fit_Data/sum(Fit_Data,'all');

        % Store RT distribution:
        model_fit_NB(sub_ind,speed_ind).Fit = Fit_Data;

        % Compute AICc:
        L = sum(log(Fit_Data(Sub_T+1)));
        model_fit_NB(sub_ind,speed_ind).AICc = 2*length(nbin)-2*L +...
            (2*length(nbin)^2+2*length(nbin))/(length(Sub_T)-length(nbin)-1);
    end
end
% save('model_fit_NB.mat','model_fit_NB')