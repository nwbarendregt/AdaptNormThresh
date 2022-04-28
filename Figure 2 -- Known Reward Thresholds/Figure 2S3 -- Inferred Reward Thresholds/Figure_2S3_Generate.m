% Figure_2S3_Generate.m
% Script used to generate Figure 2 -- Supplemental Figure 3 from 
% Barendregt et al., 2022.

clear

% Load windowed threshold data:
load('IRC_Threshold_Data.mat')

% Generate changepoint-trigged-average of thresholds for high-to-low reward
% switches (Fig. 2S2 B)
figure; hold on
for i = 1:length(m_r)
    plot(-W:dt:W,mean(thresh_HL{i},1),'linewidth',15) % Finite m_r thresholds
end
plot(-W:dt:W,mean(thresh_HL_broad,1),':','linewidth',15) % Infinite m_r thresholds
plot([0 0],[1.1 1.4],'k--','linewidth',5)

% Generate changepoint-trigged-average of thresholds for low-to-high reward
% switches (Fig. 2S2 C)
figure; hold on
for i = 1:length(m_r)
    plot(-W:dt:W,mean(thresh_LH{i},1),'linewidth',15) % Finite m_r thresholds
end
plot(-W:dt:W,mean(thresh_LH_broad,1),':','linewidth',15) % Infinite m_r thresholds
plot([0 0],[1.1 1.4],'k--','linewidth',5)