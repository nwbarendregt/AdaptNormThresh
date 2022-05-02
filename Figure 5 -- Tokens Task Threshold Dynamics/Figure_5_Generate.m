% Figure_5S1_Generate.m
% Script used to generate Figure 5 Barendregt et al., 2022.

clear

% Load motif data for slow version of tokens task:
load('Tok_Threshold_Motifs_Slow.mat')

% Generate motif colormap for slow version of tokens task (Fig. 5B):
figure
imagesc([c(1) c(end)],[R_c(1) R_c(end)],type')
set(gca,'ydir','normal')
colormap([141 211 199;...
    253 180 98;...
    190 186 218;...
    251 128 114]/255)

% Load motif data for fast version of tokens task:
load('Tok_Threshold_Motifs_Fast.mat')

% Generate motif colormap for fast version of tokens task (Fig. 5C):
figure
imagesc([c(1) c(end)],[R_c(1) R_c(end)],type')
set(gca,'ydir','normal')
colormap([141 211 199;...
    253 180 98;...
    190 186 218;...
    251 128 114]/255)

% Define simulation parameters for empirical RT distributions of each
% motif:
t_d = 0.170;
N_trial = 1e4;

% Generate motif i behavior and response distribution (Fig. 5i):
figure
thresh_g = tok_Bellmans_g(Nt,t_d,R_c(1),R_i,@(t) c(end),tol);
thresh_y = log(thresh_g./(1-thresh_g));
stairs(0:Nt,thresh_y,'linewidth',15,'color','#8dd3c7')
hold on
stairs(0:Nt,-thresh_y,'linewidth',15,'color','#8dd3c7')
line([0 Nt],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 Nt])
RT = NaN(1,N_trial);
for n = 1:N_trial
    RT(n) = tok_sim_norm(Nt,thresh_g,0);
end
RT = histcounts(RT,-0.5:(Nt+0.5),'normalization','pdf');
figure
stairs(-0.5:(Nt+0.5),[RT 0],'linewidth',15,'color','#8dd3c7')
xlim([0 Nt])

% Generate motif ii behavior and response distribution (Fig. 5ii):
figure
thresh_g = tok_Bellmans_g(Nt,t_d,R_c(10),R_i,@(t) c(1),tol);
thresh_y = log(thresh_g./(1-thresh_g));
stairs(0:Nt,thresh_y,'linewidth',15,'color','#fdb462')
hold on
stairs(0:Nt,-thresh_y,'linewidth',15,'color','#fdb462')
line([0 Nt],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 Nt])
RT = NaN(1,N_trial);
for n = 1:N_trial
    RT(n) = tok_sim_norm(Nt,thresh_g,0);
end
RT = histcounts(RT,-0.5:(Nt+0.5),'normalization','pdf');
figure
stairs(-0.5:(Nt+0.5),[RT 0],'linewidth',15,'color','#fdb462')
xlim([0 Nt])

% Generate motif iii behavior and response distribution (Fig. 5iii):
figure
thresh_g = tok_Bellmans_g(Nt,t_d,R_c(10),R_i,@(t) c(160),tol);
thresh_y = log(thresh_g./(1-thresh_g));
stairs(0:Nt,thresh_y,'linewidth',15,'color','#bebada')
hold on
stairs(0:Nt,-thresh_y,'linewidth',15,'color','#bebada')
line([0 Nt],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 Nt])
RT = NaN(1,N_trial);
for n = 1:N_trial
    RT(n) = tok_sim_norm(Nt,thresh_g,0);
end
RT = histcounts(RT,-0.5:(Nt+0.5),'normalization','pdf');
figure
stairs(-0.5:(Nt+0.5),[RT 0],'linewidth',15,'color','#bebada')
xlim([0 Nt])

% Generate motif iv behavior and response distribution (Fig. 5iv):
figure
thresh_g = tok_Bellmans_g(Nt,t_d,R_c(10),R_i,@(t) c(175),tol);
thresh_y = log(thresh_g./(1-thresh_g));
stairs(0:Nt,thresh_y,'linewidth',15,'color','#fb8072')
hold on
stairs(0:Nt,-thresh_y,'linewidth',15,'color','#fb8072')
line([0 Nt],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 Nt])
RT = NaN(1,N_trial);
for n = 1:N_trial
    RT(n) = tok_sim_norm(Nt,thresh_g,0);
end
RT = histcounts(RT,-0.5:(Nt+0.5),'normalization','pdf');
figure
stairs(-0.5:(Nt+0.5),[RT 0],'linewidth',15,'color','#fb8072')
xlim([0 Nt])