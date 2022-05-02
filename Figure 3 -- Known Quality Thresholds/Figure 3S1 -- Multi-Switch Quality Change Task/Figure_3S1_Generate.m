% Figure_3S1_Generate.m
% Script used to generate Figure 3 -- Supplemental Figure 1 from 
% Barendregt et al., 2022.

clear

T = 5; dt = 0.005;  t_i = 1;
dg = 0.001;
R = 5; c = @(t) 1; N = 5;
m = NaN(1,T/dt+1);
m(1:29) = 1.0; 
m(30:98) = 5.5;
m(99:161) = 2.8;
m(162:193) = 9.8;
m(194:end) = 8;

% Simulate reward timeseries (Fig. 3 -- Supplement 1 A) and normative 
% model (Fig. 3 -- Supplement 1 B):
SC_Threshold_Schematic(T,dt,t_i,dg,m,c,R,N)