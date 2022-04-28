% Figure_2_Generate.m
% Script used to generate Figure 2 -- Supplemental Figure 1 from 
% Barendregt et al., 2022.

% Define plotting domain and Gaussian generating function:
x = linspace(-10,10,1000);
f = @(x,mu,sigma) exp(-0.5*((x-mu)/sigma).^2)/(sigma*sqrt(2*pi));

% Plot low evidence quality distributions (Fig. 2 -- Supplement 1 A):
mu = 1; sigma = 1;
figure
plot(x,f(x,mu,sigma),'b','linewidth',15)
hold on
plot(x,f(x,-mu,sigma),'r','linewidth',15)
plot([0 0],[0 0.4],'k--','linewidth',5)

% Define simulation parameters for belief trajectories:
T = 1; dt = 0.005; N_trial = 5;

% Calculate scaled SNR m:
m = 2*mu^2/(sigma^2);

% Generate belief trajectories:
y = RDMD_trial_generate(2,T,dt,0,N_trial);

% Plot low evidence quality belief trajectories (Fig. 2 -- Supplement 1 B):
figure
plot(0:dt:T,y,'linewidth',5)
hold on
plot([0 T],[0 0],'k--','linewidth',5)

% Plot high evidence quality distributions (Fig. 2 -- Supplement 1 C):
mu = 5; sigma = 1;
figure
plot(x,f(x,mu,sigma),'b','linewidth',15)
hold on
plot(x,f(x,-mu,sigma),'r','linewidth',15)
plot([0 0],[0 0.4],'k--','linewidth',5)

% Generate and plot high evidence quality belief 
% trajectories (Fig. 2 -- Supplement 1 B):
m = 2*mu^2/(sigma^2);
y = RDMD_trial_generate(2,T,dt,0,N_trial);

figure
plot(0:dt:T,y,'linewidth',5)
hold on
plot([0 T],[0 0],'k--','linewidth',5)