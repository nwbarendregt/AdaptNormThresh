% tok_Bellmans_g.m
% Function used to obtain normative thresholds in likelihood space via 
% dynamic programming for the tokens task from Barendregt et al., 2022.

function [g_theta,rho] = tok_Bellmans_g(Nt,t_d,R_c,R_i,c,tol)

% Define time and discretizations:
N = 0:Nt; t_i = 0.5+t_d*(Nt:-1:0);
dt = 0.2; t = 0:dt:(dt*Nt);

% Define overall likelihood discretization:
g = [];
for i = 0:Nt
    for j = 0:i
        g = [g prob_choice(i-j,j,Nt)];
    end
end
g = sort(unique(g)); g_i = find(g==0.5);

% Initialize secant method:
rho = 0.1; k = 1;

% Pre-allocate value function V and maximal index of value function V_I:
V = NaN(length(g),length(N)); V_I = NaN(length(g),length(N));

% Calculate value function using backward induction:
[V(1,end),V_I(1,end)] = max([g(1)*R_c+(1-g(1))*R_i-t_i(end)*rho(k) (1-g(1))*R_c+g(1)*R_i-t_i(end)*rho(k)],[],2);
[V(end,end),V_I(end,end)] = max([g(end)*R_c+(1-g(end))*R_i-t_i(end)*rho(k) (1-g(end))*R_c+g(end)*R_i-t_i(end)*rho(k)],[],2);
for j = (length(N)-1):-1:1
    
    % Define current likelihood discretization:
    gt = [];
    for i = 0:(j-1)
        gt = [gt prob_choice((j-1)-i,i,Nt)];
    end
    gt = unique(gt);

    % Compute future value assuming singular responses at likelihoods of 0
    % and 1:
    for i = 1:length(g)
        if sum(gt==g(i))~=0
            if g(i) == 0
                [V(i,j),V_I(i,j)] = max([g(1)*R_c+(1-g(1))*R_i-t_i(j)*rho(k) (1-g(1))*R_c+g(1)*R_i-t_i(j)*rho(k) V(1,j+1)-(c(t(j))+rho(k))*dt],[],2);
            elseif g(i) == 1
                [V(i,j),V_I(i,j)] = max([g(end)*R_c+(1-g(end))*R_i-t_i(j)*rho(k) (1-g(end))*R_c+g(end)*R_i-t_i(j)*rho(k) V(end,j+1)-(c(t(j))+rho(k))*dt],[],2);
            else
                f = find(~isnan(V(:,j+1)));
                ip = f(f>=i); im = f(f<=i); ip = ip(1); im = im(end);
                [V(i,j),V_I(i,j)] = max([g(i)*R_c+(1-g(i))*R_i-t_i(j)*rho(k) (1-g(i))*R_c+g(i)*R_i-t_i(j)*rho(k) 0.5*(V(ip,j+1)+V(im,j+1))-(c(t(j))+rho(k))*dt],[],2);
            end
        end
    end
end

% Store initial value to measure convergence and perform second
% initialization of secant method:
V_rho(k) = V(g_i,1); rho = [rho 0.9]; k = 2;

V = NaN(length(g),length(N)); V_I = NaN(length(g),length(N));
[V(1,end),V_I(1,end)] = max([g(1)*R_c+(1-g(1))*R_i-t_i(end)*rho(k) (1-g(1))*R_c+g(1)*R_i-t_i(end)*rho(k)],[],2);
[V(end,end),V_I(end,end)] = max([g(end)*R_c+(1-g(end))*R_i-t_i(end)*rho(k) (1-g(end))*R_c+g(end)*R_i-t_i(end)*rho(k)],[],2);
for j = (length(N)-1):-1:1
    gt = [];
    for i = 0:(j-1)
        gt = [gt prob_choice((j-1)-i,i,Nt)];
    end
    gt = unique(gt);

    for i = 1:length(g)
        if sum(gt==g(i))~=0
            if g(i) == 0
                [V(i,j),V_I(i,j)] = max([g(1)*R_c+(1-g(1))*R_i-t_i(j)*rho(k) (1-g(1))*R_c+g(1)*R_i-t_i(j)*rho(k) V(1,j+1)-(c(t(j))+rho(k))*dt],[],2);
            elseif g(i) == 1
                [V(i,j),V_I(i,j)] = max([g(end)*R_c+(1-g(end))*R_i-t_i(j)*rho(k) (1-g(end))*R_c+g(end)*R_i-t_i(j)*rho(k) V(end,j+1)-(c(t(j))+rho(k))*dt],[],2);
            else
                f = find(~isnan(V(:,j+1)));
                ip = f(f>=i); im = f(f<=i); ip = ip(1); im = im(end);
                [V(i,j),V_I(i,j)] = max([g(i)*R_c+(1-g(i))*R_i-t_i(j)*rho(k) (1-g(i))*R_c+g(i)*R_i-t_i(j)*rho(k) 0.5*(V(ip,j+1)+V(im,j+1))-(c(t(j))+rho(k))*dt],[],2);
            end
        end
    end
end

V_rho(k) = V(g_i,1);

% Continue interating using secant method until initial value has
% sufficiently converged:
while abs(V_rho(k)) > tol
    
    % Update new reward rate using secant method:
    k = k+1; rho(k) = rho(k-1)-V_rho(k-1)*(rho(k-1)-rho(k-2))/(V_rho(k-1)-V_rho(k-2));
    
    V = NaN(length(g),length(N)); V_I = NaN(length(g),length(N));
    [V(1,end),V_I(1,end)] = max([g(1)*R_c+(1-g(1))*R_i-t_i(end)*rho(k) (1-g(1))*R_c+g(1)*R_i-t_i(end)*rho(k)],[],2);
    [V(end,end),V_I(end,end)] = max([g(end)*R_c+(1-g(end))*R_i-t_i(end)*rho(k) (1-g(end))*R_c+g(end)*R_i-t_i(end)*rho(k)],[],2);
    for j = (length(N)-1):-1:1
        gt = [];
        for i = 0:(j-1)
            gt = [gt prob_choice((j-1)-i,i,Nt)];
        end
        gt = unique(gt);

        for i = 1:length(g)
            if sum(gt==g(i))~=0
                if g(i) == 0
                    [V(i,j),V_I(i,j)] = max([g(1)*R_c+(1-g(1))*R_i-t_i(j)*rho(k) (1-g(1))*R_c+g(1)*R_i-t_i(j)*rho(k) V(1,j+1)-(c(t(j))+rho(k))*dt],[],2);
                elseif g(i) == 1
                    [V(i,j),V_I(i,j)] = max([g(end)*R_c+(1-g(end))*R_i-t_i(j)*rho(k) (1-g(end))*R_c+g(end)*R_i-t_i(j)*rho(k) V(end,j+1)-(c(t(j))+rho(k))*dt],[],2);
                else
                    f = find(~isnan(V(:,j+1)));
                    ip = f(f>=i); im = f(f<=i); ip = ip(1); im = im(end);
                    [V(i,j),V_I(i,j)] = max([g(i)*R_c+(1-g(i))*R_i-t_i(j)*rho(k) (1-g(i))*R_c+g(i)*R_i-t_i(j)*rho(k) 0.5*(V(ip,j+1)+V(im,j+1))-(c(t(j))+rho(k))*dt],[],2);
                end
            end
        end
    end
    V_rho(k) = V(g_i,1);  
end
% Return reward rate for converged thresholds:
rho = rho(end);

% Pre-allocate normative thresholds in likelihood space:
g_theta = NaN(1,length(N));

% Construct normative thresholds in likelihood space based off maximal
% index V_I:
for i =(length(g_theta)):-1:1
    if ~isempty(find(V_I(:,i)==1,1))
        if ~isempty(find(V_I(:,i)==3,1))
            g_theta(i) = g(find(V_I(:,i)==3,1,'last')+1);
        else
            g_theta(i:end) = 0.5;
        end
    else
        g_theta(i) = 1;
    end
end