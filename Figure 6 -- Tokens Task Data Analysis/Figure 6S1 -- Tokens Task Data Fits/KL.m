% KL.m
% Computes discrete KL divergence between subject data P and fitted model Q
% from Barendregt et al., 2022.

function d = KL(P,Q)

% Add small non-zero entries to fitted model:
Q(Q==0) = eps; Q = Q/sum(Q);

% Compute KL divergence:
d = P.*log(P./Q); d(isnan(d)) = 0;
d = sum(d);