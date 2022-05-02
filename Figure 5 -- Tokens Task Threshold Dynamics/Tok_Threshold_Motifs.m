% Tok_Threshold_Motifs.m
% Script used to generate data for Fig. 5B-C in Barendregt et al., 2022.

clear

% Define simulation parameters:
Nt = 15; R_i = -1; tol = 1e-5;
t_d = 0.170; % Slow version of tokens task
% t_d = 0.020; % Fast version of tokens task
c = linspace(0,1); R_c = linspace(0,10);

% Pre-allocate threshold classification matrix:
type = NaN(length(c),length(R_c));

for i = 1:length(c)
    for j = 1:length(R_c)

        % Calculate normative thresholds using Bellman's equation:
        thresh = tok_Bellmans_TL(Nt,t_d,R_c(j),R_i,@(t) c(i),tol);

        % Classify threshold motif:
        if sum(thresh==0)==length(thresh)
            type(i,j) = 1; % Motif i
        elseif sum(diff(thresh)<=0)==(length(thresh)-1)
            type(i,j) = 2; % Motif ii
        elseif sum(diff(find(sign(diff(thresh))==1))==1)==0
            type(i,j) = 3; % Motif iii
        else
            type(i,j) = 4; % Motif iv
        end
    end
end
% save('Tok_Threshold_Motifs_Slow.mat')