% RDMD_trial_generate.m
% Function that generates populations of belief trajectories for 
% continuous 2AFC tasks with static SNR given by Eq. (8) in 
% Barendregt et al., 2022.

function [y,p]= RDMD_trial_generate(m,T,dt,sigma,N_trial)

% Pre-allocate belief storage for belief with (y) and without (y_p) belief
% noise (y_p used for UGM simulations, where noise is added to filter 
% rather than input):
y = zeros(N_trial,round(T/dt)+1); y_p = zeros(N_trial,round(T/dt)+1); 
for i = 1:N_trial
    for j = 2:(T/dt+1)
        
        % Calculate white noise:
        dW = sqrt(dt)*randn;
        
        % Update observer belief:
        y(i,j) = y(i,j-1)+m*dt+sqrt(2*m)*dW+sigma*randn;
        y_p(i,j) = y(i,j-1)+m*dt+sqrt(2*m)*dW;
    end
end

% Convert LLR y_p to a likelihood to use as input to the UGM:
p = exp(y_p)./(1+exp(y_p));
end