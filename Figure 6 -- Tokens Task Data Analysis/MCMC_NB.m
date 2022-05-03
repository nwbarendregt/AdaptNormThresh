% MCMC_NB.m
% Evolves parameter chains for NB model using MCMC from 
% Barendregt et al., 2022.

function [theta,L] = MCMC_NB(theta_p,L_p,Nt,t_d,R_i,tol,Sub_T,stim,S,prior)

% Draw proposed new state from Gaussian proposal distribution:
theta = mvnrnd(theta_p,S);
if (theta(1) > prior(1)) && (theta(2) > prior(2)) &&...
        (theta(3) > prior(3)) && (theta(4) > prior(4)) % Reject new state if parameters are negative.

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

    % Add small non-zero entries to compute likelihood of proposed sampled:
    Synth_Data(Synth_Data == 0) = eps; Synth_Data = Synth_Data/sum(Synth_Data,'all');

    % Compute likelihood of proposed sample:
    L = sum(log(Synth_Data(Sub_T+1)));

    % Compute Hastings ratio and determine if proposed sample is accepted:
    h = exp(L-L_p); U = rand;
    if U > h
        theta = theta_p; L = L_p;
    end
else
    theta = theta_p; L = L_p;
end
end