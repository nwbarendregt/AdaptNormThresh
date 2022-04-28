% Figure_2S2_Generate.m
% Script used to generate Figure 2 -- Supplemental Figure 2 from 
% Barendregt et al., 2022.

clear

T = 5; dt = 0.005;  t_i = 1;
dg = 0.001;
m = 2; c = @(t) 1; N = 5;
R = NaN(1,T/dt+1); 
R(1:29) = 4.2; 
R(30:98) = 9.2;
R(99:161) = 8;
R(162:193) = 9.6;
R(194:end) = 6.6;

% Simulate reward timeseries (Fig. 2 -- Supplement 2 A) and normative 
% model (Fig. 2 -- Supplement 2 B):
RC_Threshold_Schematic(T,dt,t_i,dg,m,c,R,N)