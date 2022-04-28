% RC_Threshold_Schematic.m
% Function used to simulate normative model for reward change task from
% Barendregt et al., 2022.

function RC_Threshold_Schematic(T,dt,t_i,dg,m,c,R,N)

% Find normative thresholds using dynamic programming:
[theta,~] = RC_Bellmans(T,dt,t_i,dg,m,c,R); 
theta = theta(1:(T/5/dt+1)); % Truncate thresholds to avoid numerical artifacts.

% Pre-allocate ideal observer belief trajectories:
y = NaN(N,T/dt+1);
% Simulate ideal observer belief trajectories given by Eq. (8):
for i = 1:N
    y(i,1) = 0; k = 1;
    while abs(y(i,k)) < theta(k)
        k = k+1;
        y(i,k) = y(i,k-1)+m*dt+sqrt(2*m*dt)*randn;
    end
    % Terminate evidence accumulation when belief crosses threshold:
    y(i,k) = sign(y(i,k))*theta(k);
end

% Plot reward time series:
figure
plot(0:dt:(T/5),R(1:(T/5/dt+1)),'linewidth',15,'color','k')
% Plot belief trajectories and normative thresholds:
figure
plot(0:dt:(T/5),y(:,1:(T/5/dt+1)),'linewidth',5)
hold on
plot(0:dt:(T/5),theta,'linewidth',15,'color','k')
plot(0:dt:(T/5),-theta,'linewidth',15,'color','k')
line([0 T/5],[0 0],'linestyle','--','color','k','linewidth',5)
xlim([0 T/5])