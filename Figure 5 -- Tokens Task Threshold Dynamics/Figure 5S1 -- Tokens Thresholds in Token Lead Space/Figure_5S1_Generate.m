% Figure_5S1_Generate.m
% Script used to generate Figure 5 -- Supplemental Figure 1 from 
% Barendregt et al., 2022.

clear

% Define simulation parameters:
Nt = 15; R_i = -1; tol = 1e-5;
t_d = 0.170; % Slow version of tokens task
c = linspace(0,1,1000); R_c = linspace(0,10,1000);

% Define simulation parameters for empirical RT distributions of each
% motif:
N_trial = 1e4;

% Generate motif i behavior and response distribution (Fig. 5S1 i):
figure
thresh_TL = tok_Bellmans_TL(Nt,t_d,R_c(1),R_i,@(t) c(end),tol);
stairs(0:Nt,thresh_TL,'linewidth',15,'color','#8dd3c7')
hold on
stairs(0:Nt,-thresh_TL,'linewidth',15,'color','#8dd3c7')
line([0 Nt],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 Nt])

% Generate motif ii behavior and response distribution (Fig. 5S1 ii):
figure
thresh_TL = tok_Bellmans_TL(Nt,t_d,R_c(10),R_i,@(t) c(1),tol);
stairs(0:Nt,thresh_TL,'linewidth',15,'color','#fdb462')
hold on
stairs(0:Nt,-thresh_TL,'linewidth',15,'color','#fdb462')
line([0 Nt],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 Nt])

% Generate motif iii behavior and response distribution (Fig. 5S1 iii):
figure
thresh_TL = tok_Bellmans_TL(Nt,t_d,R_c(10),R_i,@(t) c(160),tol);
stairs(0:Nt,thresh_TL,'linewidth',15,'color','#bebada')
hold on
stairs(0:Nt,-thresh_TL,'linewidth',15,'color','#bebada')
line([0 Nt],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 Nt])

% Generate motif iv behavior and response distribution (Fig. 5S1 iv):
figure
thresh_TL = tok_Bellmans_TL(Nt,t_d,R_c(10),R_i,@(t) c(175),tol);
stairs(0:Nt,thresh_TL,'linewidth',15,'color','#fb8072')
hold on
stairs(0:Nt,-thresh_TL,'linewidth',15,'color','#fb8072')
line([0 Nt],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 Nt])