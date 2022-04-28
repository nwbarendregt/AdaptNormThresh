% RC_Bellmans.m
% Function used to obtain normative thresholds via dynamic programming for 
% reward change task from Barendregt et al., 2022.

function [theta,rho] = RC_Bellmans(T,dt,t_i,dg,m,c,R)

% Define time and likelihood discretizations:
t = 0:dt:T;
g = dg:dg:(1-dg); g_i = find(g==0.5);

% Calculate Gaussian mean from SNR m (assuming standard deviation of 1):
mu = sqrt(m/2);

% Construct likelihood transfer function given by Eq. (14):
P_gg = diag(g)*exp(-0.5/dt*(log(((g-1)'*g)./(g'*(g-1)))/(2*mu)-mu*dt).^2)/sqrt(2*pi*dt)./(2*g*mu-2*g.^2*mu)+...
    diag(1-g)*exp(-0.5/dt*(log(((g-1)'*g)./(g'*(g-1)))/(2*mu)+mu*dt).^2)/sqrt(2*pi*dt)./(2*g*mu-2*g.^2*mu);
P_gg = P_gg*diag(1./sum(P_gg,2));

% Initialize secant method:
rho = 0.1; tol = 1e-5; k = 1;

% Pre-allocate value function V and maximal index of value function V_I:
V = NaN(length(g),length(t)); V_I = NaN(length(g),length(t));

% Calculate value function using backward induction:
[V(:,end),V_I(:,end)] = max([g'*R(end)-t_i*rho(k) (1-g)'*R(end)-t_i*rho(k)],[],2);
for j = (length(t)-1):-1:1
    [V(:,j),V_I(:,j)] = max([g'*R(j)-t_i*rho(k) (1-g)'*R(j)-t_i*rho(k) P_gg*V(:,j+1)-c(t(j))*dt-rho(k)*dt],[],2);
end

% Store initial value to measure convergence and perform second
% initialization of secant method:
V_rho(k) = V(g_i,1); rho = [rho 0.9]; k = 2;

V = NaN(length(g),length(t)); V_I = NaN(length(g),length(t));
[V(:,end),V_I(:,end)] = max([g'*R(end)-t_i*rho(k) (1-g)'*R(end)-t_i*rho(k)],[],2);
for j = (length(t)-1):-1:1
    [V(:,j),V_I(:,j)] = max([g'*R(j)-t_i*rho(k) (1-g)'*R(j)-t_i*rho(k) P_gg*V(:,j+1)-c(t(j))*dt-rho(k)*dt],[],2);
end
V_rho(k) = V(g_i,1);

% Continue interating using secant method until initial value has
% sufficiently converged:
while abs(V_rho(k)) > tol

    % Update new reward rate using secant method:
    k = k+1; rho(k) = rho(k-1)-V_rho(k-1)*(rho(k-1)-rho(k-2))/(V_rho(k-1)-V_rho(k-2));
    
    V = NaN(length(g),length(t)); V_I = NaN(length(g),length(t));
    [V(:,end),V_I(:,end)] = max([g'*R(end)-t_i*rho(k) (1-g)'*R(end)-t_i*rho(k)],[],2);
    for j = (length(t)-1):-1:1
        [V(:,j),V_I(:,j)] = max([g'*R(j)-t_i*rho(k) (1-g)'*R(j)-t_i*rho(k) P_gg*V(:,j+1)-c(t(j))*dt-rho(k)*dt],[],2);
    end
    V_rho(k) = V(g_i,1); 
end
% Return reward rate for converged thresholds:
rho = rho(end);

% Pre-allocate normative thresholds in likelihood space:
g_theta = NaN(1,length(t));

% Construct normative thresholds in likelihood space based off maximal
% index V_I:
for i = length(t):-1:1
    if ~isempty(find(V_I(:,i)==1,1))
        if sum(V_I(:,i)==1)==length(g)
            g_theta(i) = 0.5;
        else
            g_theta(i) = g(find(V_I(:,i)==1,1));
        end
    else
        g_theta(i) = 1;
    end
    if g_theta(i) == 0.5
        g_theta((i+1):end) = 0.5;
    end
end

% Convert normative thresholds to LLR space:
theta = log(g_theta./(1-g_theta));