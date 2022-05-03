function p = prob_choice(nTokensPro, nTokensCon, nTokensTot)

% PROB_CHOICE Calculate the probability that a given bin will receive more tokens
%
% p = prob_choice(nTokensPro, nTokensCon, nTokensTot)
%
% nTokensPro      - Number of tokens in favor of the choice
% nTokensCon      - Number of tokens against the choice
% nTokensTot      - Total number of tokens
%
% p               - Probability that Pro will be larger than Con when all tokens are committed

if rem(nTokensTot,2)==0
   nTokensTot = nTokensTot-1;
end

nNeedPro = ceil(nTokensTot / 2) - nTokensCon;  % Need this many more tokens to be in favor

nTokensLeft = nTokensTot - nTokensPro - nTokensCon;

p = 0;

for n=0:nNeedPro-1
   
   if nTokensLeft<n
      return;
   end
   
   % pNK = the probability of getting exactly K pro given that we have N left
   nCombs = 2^nTokensLeft;
   nGoodCombs = factorial(nTokensLeft) / (factorial(n) * factorial(nTokensLeft-n));
   
   pNK = nGoodCombs / nCombs;
   
   p = p + pNK;
   
end


