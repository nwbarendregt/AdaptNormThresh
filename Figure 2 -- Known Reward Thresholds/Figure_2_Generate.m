% Figure_2_Generate.m
% Script used to generate Figure 2 from Barendregt et al., 2022.

clear

% Define simulation variables for reward change task simulations:
T = 5; dt = 0.005;  t_i = 1;
dg = 0.001;
m = 2; c = @(t) 1; N = 5;
R = NaN(1,T/dt+1); 
R(1:100) = 4; R(101:end) = 5;

% Simulate reward timeseries (Fig. 2A) and normative model (Fig. 2B):
RC_Threshold_Schematic(T,dt,t_i,dg,m,c,R,N)

% Load motif data:
load('RC_Threshold_Motifs.mat')

% Generate colormap of motifs (Fig. 2C):
figure
imagesc([R_1(1) R_1(end)],[R_2(1) R_2(end)],type')
set(gca,'ydir','normal')
colormap([253 180 98;...
    190 186 218;...
    251 128 114;...
    128 177 211]/255)

% Define simulation parameters to construct empirical response time
% distributions:
N_trial = 1e4;
mid = linspace(0,1,25);
edges = [mid-0.5*(mid(1)+mid(2)) mid(end)+0.5*(mid(1)+mid(2))];

% Generate motif i behavior and response distribution (Fig. 2i):
figure
plot(0:dt:(T/5),squeeze(thresh(1,100,:)),'linewidth',15,'color','#bebada')
hold on
plot(0:dt:(T/5),squeeze(-thresh(1,100,:)),'linewidth',15,'color','#bebada')
line([0 T/5],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 T/5])
NB_RT = NaN(1,N_trial);
y = RDMD_trial_generate(m,T/5,dt,0,N_trial);
for n = 1:N_trial
    NB_RT(n) = RDMD_sim_norm(y(n,:),T/5,dt,squeeze(thresh(1,100,:)),0);
end
NB_RT = histcounts(NB_RT,edges,'normalization','pdf');
figure
stairs(edges,[NB_RT 0],'linewidth',15,'color','#bebada')
xlim([0 T/5])

% Generate motif ii behavior and response distribution (Fig. 2ii):
figure
plot(0:dt:(T/5),squeeze(thresh(40,60,:)),'linewidth',15,'color','#fb8072')
hold on
plot(0:dt:(T/5),squeeze(-thresh(40,60,:)),'linewidth',15,'color','#fb8072')
line([0 T/5],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 T/5])
NB_RT = NaN(1,N_trial);
y = RDMD_trial_generate(m,T/5,dt,0,N_trial);
for n = 1:N_trial
    NB_RT(n) = RDMD_sim_norm(y(n,:),T/5,dt,squeeze(thresh(40,60,:)),0);
end
NB_RT = histcounts(NB_RT,edges,'normalization','pdf');
figure
stairs(edges,[NB_RT 0],'linewidth',15,'color','#fb8072')
xlim([0 T/5])

% Generate motif iii behavior and response distribution (Fig. 2iii):
figure
plot(0:dt:(T/5),squeeze(thresh(50,50,:)),'linewidth',15,'color','#80b1d3')
hold on
plot(0:dt:(T/5),squeeze(-thresh(50,50,:)),'linewidth',15,'color','#80b1d3')
line([0 T/5],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 T/5])
NB_RT = NaN(1,N_trial);
y = RDMD_trial_generate(m,T/5,dt,0,N_trial);
for n = 1:N_trial
    NB_RT(n) = RDMD_sim_norm(y(n,:),T/5,dt,squeeze(thresh(50,50,:)),0);
end
NB_RT = histcounts(NB_RT,edges,'normalization','pdf');
figure
stairs(edges,[NB_RT 0],'linewidth',15,'color','#80b1d3')
xlim([0 T/5])

% Generate motif iv behavior and response distribution (Fig. 2iv):
figure
plot(0:dt:(T/5),squeeze(thresh(55,45,:)),'linewidth',15,'color','#fb8072')
hold on
plot(0:dt:(T/5),squeeze(-thresh(55,45,:)),'linewidth',15,'color','#fb8072')
line([0 T/5],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 T/5])
NB_RT = NaN(1,N_trial);
y = RDMD_trial_generate(m,T/5,dt,0,N_trial);
for n = 1:N_trial
    NB_RT(n) = RDMD_sim_norm(y(n,:),T/5,dt,squeeze(thresh(55,45,:)),0);
end
NB_RT = histcounts(NB_RT,edges,'normalization','pdf');
figure
stairs(edges,[NB_RT 0],'linewidth',15,'color','#fb8072')
xlim([0 T/5])

% Generate motif v behavior and response distribution (Fig. 2v):
figure
plot(0:dt:(T/5),squeeze(thresh(100,1,:)),'linewidth',15,'color','#fdb462')
hold on
plot(0:dt:(T/5),squeeze(-thresh(100,1,:)),'linewidth',15,'color','#fdb462')
line([0 T/5],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 T/5])
NB_RT = NaN(1,N_trial);
y = RDMD_trial_generate(m,T/5,dt,0,N_trial);
for n = 1:N_trial
    NB_RT(n) = RDMD_sim_norm(y(n,:),T/5,dt,squeeze(thresh(100,1,:)),0);
end
NB_RT = histcounts(NB_RT,edges,'normalization','pdf');
figure
stairs(edges,[NB_RT 0],'linewidth',15,'color','#fdb462')
xlim([0 T/5])