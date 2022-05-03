% Figure_6S1_Generate.m
% Script used to generate Figure 6 -- Supplemental Figure 1 
% from Barendregt et al., 2022.

clear

% Load subject and MCMC data:
load('trials.mat'); subs = 20; speeds = 2;
load('model_fit_NB.mat');
load('model_fit_Const.mat');
load('model_fit_UGM.mat');

% Define parameters to simulate models:
Nt = 15;

% Format data from NB model:
NB = model_fit_NB(:);
model_AICc = NaN(1,length(NB));

% Find best-fit and worst-fit NB model:
for i = 1:length(NB)
    model_AICc(i) = NB(i).AICc;
end
[~,m] = min(model_AICc); [~,M] = max(model_AICc);

% Plot best-fit model and associated subject data:
Sub_ID = NB(m).idSubject;
speed = NB(m).speed; t_d = speed/1000;
Sub_T = trials.nDecisionToken((trials.nSpeedFast == speed) & (trials.idSubject == Sub_ID));
Sub_Data = histcounts(Sub_T,-0.5:1:(Nt+0.5),'normalization','probability');
figure
stairs(-0.5:Nt+0.5,[Sub_Data 0],'k','linewidth',15)
hold on; xlim([0 Nt])
stairs(-0.5:Nt+0.5,[NB(m).Fit 0],'linewidth',15,'color','#1b9e77');

% Display KL divergence between subject data and fitted model:
disp(KL(Sub_Data,NB(m).Fit))

% Plot worst-fit model and associated subject data:
Sub_ID = NB(M).idSubject;
speed = NB(M).speed; t_d = speed/1000;
Sub_T = trials.nDecisionToken((trials.nSpeedFast == speed) & (trials.idSubject == Sub_ID));
Sub_Data = histcounts(Sub_T,-0.5:1:(Nt+0.5),'normalization','probability');
figure
stairs(-0.5:Nt+0.5,[Sub_Data 0],'k','linewidth',15)
hold on; xlim([0 Nt])
stairs(-0.5:Nt+0.5,[NB(M).Fit 0],'linewidth',15,'color','#1b9e77');

% Display KL divergence between subject data and fitted model:
disp(KL(Sub_Data,NB(M).Fit))

% Format data from Const model:
Const = model_fit_Const(:);
model_AICc = NaN(1,length(Const));

% Find best-fit and worst-fit Const model:
for i = 1:length(Const)
    model_AICc(i) = Const(i).AICc;
end
[~,m] = min(model_AICc); [~,M] = max(model_AICc);

% Plot best-fit model and associated subject data:
Sub_ID = Const(m).idSubject;
speed = Const(m).speed; t_d = speed/1000;
Sub_T = trials.nDecisionToken((trials.nSpeedFast == speed) & (trials.idSubject == Sub_ID));
Sub_Data = histcounts(Sub_T,-0.5:1:(Nt+0.5),'normalization','probability');
figure
stairs(-0.5:Nt+0.5,[Sub_Data 0],'k','linewidth',15)
hold on; xlim([0 Nt])
stairs(-0.5:Nt+0.5,[Const(m).Fit 0],'linewidth',15,'color','#d95f02');

% Display KL divergence between subject data and fitted model:
disp(KL(Sub_Data,Const(m).Fit))

% Plot worst-fit model and associated subject data:
Sub_ID = Const(M).idSubject;
speed = Const(M).speed; t_d = speed/1000;
Sub_T = trials.nDecisionToken((trials.nSpeedFast == speed) & (trials.idSubject == Sub_ID));
Sub_Data = histcounts(Sub_T,-0.5:1:(Nt+0.5),'normalization','probability');
figure
stairs(-0.5:Nt+0.5,[Sub_Data 0],'k','linewidth',15)
hold on; xlim([0 Nt])
stairs(-0.5:Nt+0.5,[Const(M).Fit 0],'linewidth',15,'color','#d95f02');

% Display KL divergence between subject data and fitted model:
disp(KL(Sub_Data,Const(M).Fit))

% Format data from UGM:
UGM = model_fit_UGM(:);
model_AICc = NaN(1,length(UGM));

% Find best-fit and worst-fit UGM:
for i = 1:length(UGM)
    model_AICc(i) = UGM(i).AICc;
end
[~,m] = min(model_AICc); [~,M] = max(model_AICc);

% Plot best-fit model and associated subject data:
Sub_ID = UGM(m).idSubject;
speed = UGM(m).speed; t_d = speed/1000;
Sub_T = trials.nDecisionToken((trials.nSpeedFast == speed) & (trials.idSubject == Sub_ID));
Sub_Data = histcounts(Sub_T,-0.5:1:(Nt+0.5),'normalization','probability');
figure
stairs(-0.5:Nt+0.5,[Sub_Data 0],'k','linewidth',15)
hold on; xlim([0 Nt])
stairs(-0.5:Nt+0.5,[UGM(m).Fit 0],'linewidth',15,'color','#7570b3');

% Display KL divergence between subject data and fitted model:
disp(KL(Sub_Data,UGM(m).Fit))

% Plot worst-fit model and associated subject data:
Sub_ID = UGM(M).idSubject;
speed = UGM(M).speed; t_d = speed/1000;
Sub_T = trials.nDecisionToken((trials.nSpeedFast == speed) & (trials.idSubject == Sub_ID));
Sub_Data = histcounts(Sub_T,-0.5:1:(Nt+0.5),'normalization','probability');
figure
stairs(-0.5:Nt+0.5,[Sub_Data 0],'k','linewidth',15)
hold on; xlim([0 Nt])
stairs(-0.5:Nt+0.5,[UGM(M).Fit 0],'linewidth',15,'color','#7570b3');

% Display KL divergence between subject data and fitted model:
disp(KL(Sub_Data,UGM(M).Fit))