% Figure_3_Generate.m
% Script used to generate Figure 3 from Barendregt et al., 2022.

clear

% Define simulation variables for SNR change task simulations:
T = 5; dt = 0.005; t_i = 1;
dg = 0.001;
R = 5; c = @(t) 1; N = 5;
m = NaN(1,T/dt+1);
m(1:100) = 5; m(101:end) = 2;

% Simulate reward timeseries (Fig. 3A) and normative model (Fig. 3B):
SC_Threshold_Schematic(T,dt,t_i,dg,m,c,R,N)

% Load motif data:
load('SC_Threshold_Motifs.mat')

% Generate colormap of motifs (Fig. 3C):
figure
imagesc([m_1(1) m_1(end)],[m_2(1) m_2(end)],type')
set(gca,'ydir','normal')
colormap([128 177 211;...
    253 180 98;...
    255 237 111]/255)

% Define simulation parameters to construct empirical response time
% distributions:
N_trial = 1e4;
mid = linspace(0,1,25);
edges = [mid-0.5*(mid(1)+mid(2)) mid(end)+0.5*(mid(1)+mid(2))];

% Generate motif i behavior and response distribution (Fig. 3i):
figure
plot(0:dt:(T/5),squeeze(thresh(1,100,:)),'linewidth',15,'color','#ffed6f')
hold on
plot(0:dt:(T/5),squeeze(-thresh(1,100,:)),'linewidth',15,'color','#ffed6f')
line([0 T/5],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 T/5])
NB_RT = NaN(1,N_trial);
m = NaN(1,T/dt+1); m(1:100) = m_1(1); m(101:end) = m_2(100);
y = SC_RDMD_trial_generate(m,T/5,dt,0,N_trial);
for n = 1:N_trial
    NB_RT(n) = RDMD_sim_norm(y(n,:),T/5,dt,squeeze(thresh(1,100,:)),0);
end
NB_RT = histcounts(NB_RT,edges,'normalization','pdf');
figure
stairs(edges,[NB_RT 0],'linewidth',15,'color','#ffed6f')
xlim([0 T/5])

% Generate motif ii behavior and response distribution (Fig. 3ii):
figure
plot(0:dt:(T/5),squeeze(thresh(50,50,:)),'linewidth',15,'color','#80b1d3')
hold on
plot(0:dt:(T/5),squeeze(-thresh(50,50,:)),'linewidth',15,'color','#80b1d3')
line([0 T/5],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 T/5])
NB_RT = NaN(1,N_trial);
m = NaN(1,T/dt+1); m(1:100) = m_1(50); m(101:end) = m_2(50);
y = SC_RDMD_trial_generate(m,T/5,dt,0,N_trial);
for n = 1:N_trial
    NB_RT(n) = RDMD_sim_norm(y(n,:),T/5,dt,squeeze(thresh(50,50,:)),0);
end
NB_RT = histcounts(NB_RT,edges,'normalization','pdf');
figure
stairs(edges,[NB_RT 0],'linewidth',15,'color','#80b1d3')
xlim([0 T/5])

% Generate motif iii behavior and response distribution (Fig. 3iii):
figure
plot(0:dt:(T/5),squeeze(thresh(100,1,:)),'linewidth',15,'color','#fdb462')
hold on
plot(0:dt:(T/5),squeeze(-thresh(100,1,:)),'linewidth',15,'color','#fdb462')
line([0 T/5],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 T/5])
NB_RT = NaN(1,N_trial);
m = NaN(1,T/dt+1); m(1:100) = m_1(100); m(101:end) = m_2(1);
y = SC_RDMD_trial_generate(m,T/5,dt,0,N_trial);
for n = 1:N_trial
    NB_RT(n) = RDMD_sim_norm(y(n,:),T/5,dt,squeeze(thresh(100,1,:)),0);
end
NB_RT = histcounts(NB_RT,edges,'normalization','pdf');
figure
stairs(edges,[NB_RT 0],'linewidth',15,'color','#fdb462')
xlim([0 T/5])