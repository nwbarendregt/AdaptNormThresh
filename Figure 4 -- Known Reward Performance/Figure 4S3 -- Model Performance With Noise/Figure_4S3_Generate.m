% Figure_4S3_Generate.m
% Script used to generate Figure 4 -- Supplemental Figure 3 from 
% Barendregt et al., 2022.

% Load model reward rate data from noise sweep:
load('RC_Noise_Sweep_Data.mat')

% Plot model reward rates for low-to-high reward switch (Fig. 4S3 A,B)
figure
imagesc([sigma(1) sigma(end)],[mn(1) mn(end)],NB_RR(:,:,1)'); colorbar
set(gca,'ydir','normal')
figure
imagesc([sigma(1) sigma(end)],[mn(1) mn(end)],Const_RR(:,:,1)'); colorbar
set(gca,'ydir','normal')
figure
imagesc([sigma(1) sigma(end)],[mn(1) mn(end)],UGM_RR(:,:,1)'); colorbar
set(gca,'ydir','normal')
figure
plot(sigma,NB_RR(:,50,1),'linewidth',15,'color','#1b9e77')
hold on
plot(sigma,Const_RR(:,50,1),'linewidth',15,'color','#d95f02')
plot(sigma,UGM_RR(:,50,1),'linewidth',15,'color','#7570b3')

% Plot model reward rates for low-to-high reward switch (Fig. 4S3 C,D)
figure
imagesc([sigma(1) sigma(end)],[mn(1) mn(end)],NB_RR(:,:,2)'); colorbar
set(gca,'ydir','normal')
figure
imagesc([sigma(1) sigma(end)],[mn(1) mn(end)],Const_RR(:,:,2)'); colorbar
set(gca,'ydir','normal')
figure
imagesc([sigma(1) sigma(end)],[mn(1) mn(end)],UGM_RR(:,:,2)'); colorbar
set(gca,'ydir','normal')
figure
plot(sigma,NB_RR(:,50,2),'linewidth',15,'color','#1b9e77')
hold on
plot(sigma,Const_RR(:,50,2),'linewidth',15,'color','#d95f02')
plot(sigma,UGM_RR(:,50,2),'linewidth',15,'color','#7570b3')