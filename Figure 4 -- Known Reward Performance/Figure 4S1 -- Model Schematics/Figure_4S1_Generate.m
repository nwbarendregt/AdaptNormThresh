% Figure_4S1_Generate.m
% Script used to generate Figure 4 -- Supplemental Figure 1 from 
% Barendregt et al., 2022.

clear

% Define simulation variables for reward change task simulations:
T = 1; dt = 0.005; t = 0:dt:T;
m = 10;

% Define simulation parameters to construct empirical response time
% distributions:
N_trial = 1e4;
mid = linspace(0,1,25);
edges = [mid-0.5*(mid(1)+mid(2)) mid(end)+0.5*(mid(1)+mid(2))];

% Define Const model parameters and plotting colors:
theta_0 = linspace(1,3,4);
color = {'#8d3e01','#c05402','#f26a02','#fd862a'};

f1 = figure; f2 = figure;
for i = 1:length(theta_0)

    % Construct threshold for Const model:
    thresh = theta_0(i)*ones(1,length(t));
    
    % Plot threshold dynamics of Const model (Fig. 4S1 A):
    figure(f1)
    plot(t,thresh,'linewidth',15,'color',color{i})
    hold on
    plot(t,-thresh,'linewidth',15,'color',color{i})
    ylim([-3 3])
    
    % Generate response distribution for Const model with given threshold 
    % (Fig. 4S1 A):
    RT = NaN(1,N_trial);
    y = RDMD_trial_generate(m,T,dt,0,N_trial);
    for n = 1:N_trial
        RT(n) = RDMD_sim_norm(y(n,:),T,dt,thresh,0);
    end
    RT = histcounts(RT,edges,'normalization','pdf');
    
    figure(f2)
    stairs(edges,[RT 0],'linewidth',15,'color',color{i})
    hold on
    xlim([0 T])
end
figure(f1)
plot(t,zeros(1,length(t)),'k--','linewidth',5)

% Define UGM parameters and plotting colors:
tau = 0.2; sigma = 0.7; a = 3;
theta_0 = linspace(0.1,2,4);
color = {'#4f4a8c','#655faa','#8581bc','#a5a2ce'};

f3 = figure; f4 = figure;
for i = 1:length(theta_0)

    % Construct threshold for UGM:
    thresh = theta_0(i)./(a*t);
    
    % Plot threshold dynamics of UGM (Fig. 4S1 B):
    figure(f3)
    plot(t,thresh,'linewidth',15,'color',color{i})
    hold on
    plot(t,-thresh,'linewidth',15,'color',color{i})
    xlim([0 T]); ylim([-3 3])
    
    % Generate response distribution for Const model with given threshold 
    % (Fig. 4S1 B):
    RT = NaN(1,N_trial);
    [~,p] = RDMD_trial_generate(m,T,dt,0,N_trial);
    for n = 1:N_trial
        RT(n) = RDMD_sim_UGM(p(n,:),T,dt,theta_0(i),a,sigma,tau,0);
    end
    RT = histcounts(RT,edges,'normalization','pdf');
    
    figure(f4)
    stairs(edges,[RT 0],'linewidth',15,'color',color{i})
    hold on
    xlim([0 T])
end
figure(f3)
plot(t,zeros(1,length(t)),'k--','linewidth',5)