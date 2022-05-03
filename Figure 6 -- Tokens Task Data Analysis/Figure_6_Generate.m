% Figure_6_Generate.m
% Script used to generate Figure 6 Barendregt et al., 2022.

clear

% Load MCMC and subject data:
load('trials.mat'); subs = 20; speeds = 2;
load('model_fit_NB.mat');
load('model_fit_Const.mat');
load('model_fit_UGM.mat');

% Define parameters for simulating models:
Nt = 15; R_i = -1; tol = 1e-5;

% Perform model selection using AICc and RMSE (Fig. 6A,D):
for speed_ind = 1:speeds
    model_class = zeros(3,2);
    NB_err = NaN(subs,1); Const_err = NaN(subs,1); UGM_err = NaN(subs,1);
    for sub_ind = 1:subs

        % Compare model AICc:
        model_AICc = [model_fit_NB(sub_ind,speed_ind).AICc model_fit_Const(sub_ind,speed_ind).AICc ...
            model_fit_UGM(sub_ind,speed_ind).AICc];
        [~,I] = min(model_AICc);
        model_class(I,1) = model_class(I,1)+1;
        
        % Load and format stimulus data:
        Sub_ID = model_fit_NB(sub_ind,speed_ind).idSubject;
        speed = model_fit_NB(sub_ind,speed_ind).speed; t_d = speed/1000;
        Sub_T = trials.nDecisionToken((trials.nSpeedFast == speed) & (trials.idSubject == Sub_ID));
        Sub_stim = trials.sTokenDirs((trials.nSpeedFast == speed) & (trials.idSubject == Sub_ID));
        stim = NaN(length(Sub_stim),Nt);
        for i = 1:length(Sub_stim)
            stim(i,:) = str2num(strtrim(regexprep(Sub_stim{i},'.{1}','$0 ')));
            stim(i,:) = 2*(stim(i,:)-1)-1;
        end
        
        % Compute average RMSE for NB model:
        Fit_response = NaN(50,length(Sub_T));
        MLE_NB = model_fit_NB(sub_ind,speed_ind).MLE;
        thresh_g = tok_Bellmans_g(Nt,t_d,MLE_NB(1),R_i,@(t) MLE_NB(2),tol);
        for i = 1:50
            for j = 1:length(Sub_T)
                T = tok_sim_norm(Nt,thresh_g,MLE_NB(3),stim(j,:));
                Fit_response(i,j) = round(MLE_NB(4)*randn+T);
                while (Fit_response(i,j) > Nt) || (Fit_response(i,j) < 0)
                    Fit_response(i,j) = round(MLE_NB(4)*randn+T);
                end
            end
        end
        Fit_response = mean(Fit_response,1);
        NB_err(sub_ind) = sqrt(sum((Fit_response-Sub_T).^2)/length(Sub_T));

        % Compute average RMSE for Const model:
        Fit_response = NaN(50,length(Sub_T));
        MLE_Const = model_fit_Const(sub_ind,speed_ind).MLE;
        thresh_g = MLE_Const(1)*ones(1,Nt+1);
        for i = 1:50
            for j = 1:length(Sub_T)
                T = tok_sim_norm(Nt,thresh_g,MLE_Const(2),stim(j,:));
                Fit_response(i,j) = round(MLE_Const(3)*randn+T);
                while (Fit_response(i,j) > Nt) || (Fit_response(i,j) < 0)
                    Fit_response(i,j) = round(MLE_Const(3)*randn+T);
                end
            end
        end
        Fit_response = mean(Fit_response,1);
        Const_err(sub_ind) = sqrt(sum((Fit_response-Sub_T).^2)/length(Sub_T));
        
        % Compute average RMSE for UGM:
        Fit_response = NaN(50,length(Sub_T));
        MLE_UGM = model_fit_UGM(sub_ind,speed_ind).MLE;
        for i = 1:50
            for j = 1:length(Sub_T)
                T = tok_sim_UGM(Nt,MLE_UGM(1),MLE_UGM(2),MLE_UGM(3),MLE_UGM(4),stim(j,:));
                Fit_response(i,j) = round(MLE_UGM(5)*randn+T);
                while (Fit_response(i,j) > Nt) || (Fit_response(i,j) < 0)
                    Fit_response(i,j) = round(MLE_UGM(5)*randn+T);
                end
            end
        end
        Fit_response = mean(Fit_response,1);
        UGM_err(sub_ind) = sqrt(sum((Fit_response-Sub_T).^2)/length(Sub_T));
        
        % Compare model average RMSE:
        model_err = [NB_err(sub_ind) Const_err(sub_ind) UGM_err(sub_ind)];
        [~,I] = min(model_err);
        model_class(I,2) = model_class(I,2)+1;
    end

    % Plot model selection results:
    figure
    bar(model_class')
end

% Compare subject mean RT to predicted RT from models (Fig. 6B,E):
for speed_ind = 1:speeds
    Sub_RT = NaN(subs,1); NB_RT = NaN(subs,1); Const_RT = NaN(subs,1); UGM_RT = NaN(subs,1);
    model_class = zeros(3,subs);
    for sub_ind = 1:subs
        
        % Calculate subject mean RT:
        Sub_ID = model_fit_NB(sub_ind,speed_ind).idSubject;
        speed = model_fit_NB(sub_ind,speed_ind).speed; t_d = speed/1000;
        Sub_T = trials.nDecisionToken((trials.nSpeedFast == speed) & (trials.idSubject == Sub_ID));
        Sub_RT(sub_ind) = mean(Sub_T);
        
        % Calculate model mean RT:
        NB_RT(sub_ind) = sum((0:Nt).*model_fit_NB(sub_ind,speed_ind).Fit);
        Const_RT(sub_ind) = sum((0:Nt).*model_fit_Const(sub_ind,speed_ind).Fit);
        UGM_RT(sub_ind) = sum((0:Nt).*model_fit_UGM(sub_ind,speed_ind).Fit);
        
        % Sort predictions by whether or not model was also selected using
        % AICc:
        model_AICc = [model_fit_NB(sub_ind,speed_ind).AICc model_fit_Const(sub_ind,speed_ind).AICc ...
            model_fit_UGM(sub_ind,speed_ind).AICc];
        [~,I] = min(model_AICc);
        model_class(I,sub_ind) = 1;
    end

    % Compute variance in models preditions from subject data:
    NB_RT_var = var(NB_RT-Sub_RT);
    Const_RT_var = var(Const_RT-Sub_RT);
    UGM_RT_var = var(UGM_RT-Sub_RT);
    
    % Plot mean RT results:
    figure
    plot(0:Nt,0:Nt,'k--','linewidth',5)
    hold on
    scatter(Sub_RT(model_class(1,:) == 0),NB_RT(model_class(1,:) == 0),...
        1500,'filled','MarkerEdgeColor','k','MarkerFaceColor','#1b9e77','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
    scatter(Sub_RT(model_class(2,:) == 0),Const_RT(model_class(2,:) == 0), ...
        1500,'filled','MarkerEdgeColor','k','MarkerFaceColor','#d95f02','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
    scatter(Sub_RT(model_class(3,:) == 0),UGM_RT(model_class(3,:) == 0), ...
        1500,'filled','MarkerEdgeColor','k','MarkerFaceColor','#7570b3','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
    scatter(Sub_RT(model_class(1,:) == 1),NB_RT(model_class(1,:) == 1),...
        1500,'filled','MarkerEdgeColor','k','MarkerFaceColor','#1b9e77')
    scatter(Sub_RT(model_class(2,:) == 1),Const_RT(model_class(2,:) == 1),...
        1500,'filled','MarkerEdgeColor','k','MarkerFaceColor','#d95f02')
    scatter(Sub_RT(model_class(3,:) == 1),UGM_RT(model_class(3,:) == 1), ...
        1500,'filled','MarkerEdgeColor','k','MarkerFaceColor','#7570b3')
    xlabel('Subject Mean RT')
    ylabel('Model Mean RT')

    % Display prediction variances:
    disp([NB_RT_var Const_RT_var UGM_RT_var])
end

% Classify motif of fitted normative model for each subject (Fig. 6C,F):
for speed_ind = 1:speeds
    c_MLE = NaN(subs,1); R_c_MLE = NaN(subs,1);
    type = NaN(1,subs);
    model_class = zeros(3,subs);
    for sub_ind = 1:subs

        % Load MLE parameters from fitted NB model:
        c_MLE(sub_ind) = model_fit_NB(sub_ind,speed_ind).MLE(2);
        R_c_MLE(sub_ind) = model_fit_NB(sub_ind,speed_ind).MLE(1);
        
        % Construct fitted NB threshold:
        speed = model_fit_NB(sub_ind,speed_ind).speed; t_d = speed/1000;
        thresh_TL = tok_Bellmans_TL(Nt,t_d,R_c_MLE(sub_ind),R_i,@(t) c_MLE(sub_ind),tol);

        % Classify threshold motif:
        if sum(thresh_TL==0)==length(thresh_TL)
            type(sub_ind) = 1; % Motif i
        elseif sum(diff(thresh_TL)<=0)==(length(thresh_TL)-1)
            type(sub_ind) = 2; % Motif ii
        elseif sum(diff(find(sign(diff(thresh_TL))==1))==1)==0
            type(sub_ind) = 3; % Motif iii
        else
            type(sub_ind) = 4; % Motif iv
        end

        % Sort classification by whether or not NB model was selected using
        % AICc:
        model_AICc = [model_fit_NB(sub_ind,speed_ind).AICc model_fit_Const(sub_ind,speed_ind).AICc ...
            model_fit_UGM(sub_ind,speed_ind).AICc];
        [~,I] = min(model_AICc);
        model_class(I,sub_ind) = 1;
    end

    % Plot MLEs in parameter space, color-coded by threshold motif:
    figure
    scatter(c_MLE((type == 1) & (model_class(1,:) == 1)),R_c_MLE((type == 1) & (model_class(1,:) == 1)),1500,'filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','#8dd3c7')
    hold on
    scatter(c_MLE((type == 1) & (model_class(1,:) == 0)),R_c_MLE((type == 1) & (model_class(1,:) == 1)),1500,'filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','#8dd3c7','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
    scatter(c_MLE((type == 2) & (model_class(1,:) == 1)),R_c_MLE((type == 2) & (model_class(1,:) == 1)),1500,'filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','#fdb462')
    scatter(c_MLE((type == 2) & (model_class(1,:) == 0)),R_c_MLE((type == 2) & (model_class(1,:) == 0)),1500,'filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','#fdb462','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
    scatter(c_MLE((type == 3) & (model_class(1,:) == 1)),R_c_MLE((type == 3) & (model_class(1,:) == 1)),1500,'filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','#bebada')
    scatter(c_MLE((type == 3) & (model_class(1,:) == 0)),R_c_MLE((type == 3) & (model_class(1,:) == 0)),1500,'filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','#bebada','MarkerFaceAlpha',0.25,'MarkerEdgeAlpha',0.25)
    scatter(c_MLE((type == 4) & (model_class(1,:) == 1)),R_c_MLE((type == 4) & (model_class(1,:) == 1)),1500,'filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','#fb8072')
    scatter(c_MLE((type == 4) & (model_class(1,:) == 0)),R_c_MLE((type == 4) & (model_class(1,:) == 0)),1500,'filled',...
        'MarkerEdgeColor','k','MarkerFaceColor','#fb8072')
end