% Figure_4S2_Generate.m
% Script used to generate Figure 4 -- Supplemental Figure 2 from 
% Barendregt et al., 2022.

clear

% Load noisy model performance for high-to-low reward change task:
load('RC_Fixed_Noise_Strength_R1_8.mat')

% Plot model performance for varied noise strenghts (Fig. 4S2 A):
figure
plot(noise_strength,mean(NB_RR,2),'linewidth',5,'color','#1b9e77','marker','s',...
    'markerfacecolor','auto','markersize',24)
hold on
plot(noise_strength(1),mean(NB_RR(1,:)),'marker','s','markeredgecolor','#105d46',...
    'markerfacecolor','#105d46','markersize',36)
plot(noise_strength(6),mean(NB_RR(6,:)),'marker','s','markeredgecolor','#1b9e77',...
    'markerfacecolor','#1b9e77','markersize',36)
plot(noise_strength(end),mean(NB_RR(end,:)),'marker','s','markeredgecolor','#2bdba6',...
    'markerfacecolor','#2bdba6','markersize',36)

plot(noise_strength,mean(Const_RR,2),'linewidth',5,'color','#d95f02','marker','s',...
    'markerfacecolor','auto','markersize',24)
plot(noise_strength(1),mean(Const_RR(1,:)),'marker','s','markeredgecolor','#8d3e01',...
    'markerfacecolor','#8d3e01','markersize',36)
plot(noise_strength(6),mean(Const_RR(6,:)),'marker','s','markeredgecolor','#d95f02',...
    'markerfacecolor','#d95f02','markersize',36)
plot(noise_strength(end),mean(Const_RR(end,:)),'marker','s','markeredgecolor','#fd862a',...
    'markerfacecolor','#fd862a','markersize',36)

plot(noise_strength,mean(UGM_RR,2),'linewidth',5,'color','#7570b3','marker','s',...
    'markerfacecolor','auto','markersize',24)
plot(noise_strength(1),mean(UGM_RR(1,:)),'marker','s','markeredgecolor','#4f4a8c',...
    'markerfacecolor','#4f4a8c','markersize',36)
plot(noise_strength(6),mean(UGM_RR(6,:)),'marker','s','markeredgecolor','#7570b3',...
    'markerfacecolor','#7570b3','markersize',36)
plot(noise_strength(end),mean(UGM_RR(end,:)),'marker','s','markeredgecolor','#a5a2ce',...
    'markerfacecolor','#a5a2ce','markersize',36)

% Define noise strengths for empirical RT distributions:
noise_strength = 0:0.5:1;

% Define simulation parameters for empirical RT distributions for each
% model:
f1 = figure; f2 = figure; f3 = figure;
mid = linspace(0,1,25);
edges = [mid-0.5*(mid(1)+mid(2)) mid(end)+0.5*(mid(1)+mid(2))];

for i = 1:length(noise_strength)

    % Pre-allocate RT distributions:
    NB_RT = zeros(1,length(mid)); Const_RT = zeros(1,length(mid)); UGM_RT = zeros(1,length(mid));

    % Sample strength of each noise source given fixed total noise level:
    z = noise_strength(i)*(sigma_max+mn_max);
    sigma = (min([sigma_max,z])-max([z-mn_max,0]))*rand(1,N_sub)+max([z-mn_max,0]);
    mn = z-sigma;

    for j = 1:N_sub
        RT_data = NaN(3,N_trial);

        % Generate trial data for synthetic subject:
        [y,p] = RDMD_trial_generate(m,T/5,dt,sigma(j),N_trial);
        for n = 1:N_trial

            % Calculate RT distributions for synthetic subject:
            RT_data(1,n) = RDMD_sim_norm(y(n,:),T/5,dt,NB_thresh,mn(j));
            RT_data(2,n) = RDMD_sim_norm(y(n,:),T/5,dt,Const_thresh,mn(j));
            RT_data(3,n) = RDMD_sim_UGM(p(n,:),T/5,dt,UGM_theta(1),UGM_theta(2),UGM_theta(3)+sigma(j),UGM_theta(4),mn(j));
        end

        % Store RT distribution data:
        NB_RT = NB_RT+histcounts(RT_data(1,:),edges,'normalization','pdf');
        Const_RT = Const_RT+histcounts(RT_data(2,:),edges,'normalization','pdf');
        UGM_RT = UGM_RT+histcounts(RT_data(3,:),edges,'normalization','pdf');
    end

    % Calculate RT distributions averaged over synthetic subjects:
    NB_RT = NB_RT/trapz(mid,NB_RT);
    Const_RT = Const_RT/trapz(mid,Const_RT);
    UGM_RT = UGM_RT/trapz(mid,UGM_RT);
    
    % Plot empirical RT distributions for each model (Fig. 4S2 B-D):
    figure(f1)
    stairs(edges,[NB_RT 0],'linewidth',5); xlim([0 1])
    hold on
    figure(f2)
    stairs(edges,[Const_RT 0],'linewidth',5); xlim([0 1])
    hold on
    figure(f3)
    stairs(edges,[UGM_RT 0],'linewidth',5); xlim([0 1])
    hold on
end

% Plot normative thresholds for given reward timeseries (Fig. 4S2 B inset):
figure
plot(0:dt:(T/5),squeeze(NB_thresh(1:(T/5/dt+1))),'linewidth',15,'color','#1b9e77')
hold on
plot(0:dt:(T/5),squeeze(-NB_thresh(1:(T/5/dt+1))),'linewidth',15,'color','#1b9e77')
line([0 T/5],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 T/5])