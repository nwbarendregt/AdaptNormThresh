% tok_sim_norm.m
% Generates response times for tokens task using the NB or Const model
% with added belief noise from Barendregt et al., 2022.

function T = tok_sim_norm(Nt,thresh_g,sigma,varargin)

% Use specified token movement stimulus if specified. Otherwise, generate
% random token movement stimulus:
if isempty(varargin)
    xi = [NaN 2*binornd(1,0.5,1,Nt)-1];
else
    xi = [NaN varargin{1}];
end

% Pre-allocate noise-free response time:
T = NaN;

% Convert threshold from likelihood space to LLR space:
thresh = log(thresh_g./(1-thresh_g));

% Instantly respond with random choice if threshold is trivally-collapsed:
if thresh(1) == 0
    T = 0;
else

    % Instantly respond with choice if belief is very noisy:
    g = randn*sigma;
    if g >= thresh(1)
        T = 0;
    else

        % Determine threshold crossing time and associated choice:
        k = 1;
        while isnan(T)
            k = k+1;
            g = prob_choice(sum(xi(1:k)==1),sum(xi(1:k)==-1),Nt); g = log(g./(1-g))+randn*sigma;
            if g>=thresh(k)
                T = k-1;
            elseif g<=-thresh(k)
                T = k-1;
            end
        end
    end
end