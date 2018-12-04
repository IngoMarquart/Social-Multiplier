%%%
% % BonacichCentrality 
% Calculate different Bonacich (1987) centralities
% @param: alpha - Provide alpha term
% @param: beta - Provide beta term
% @param: normAlpha - 1/0 - whether to norm Alpha term
% @param: normBeta - 1/0 - whether to set beta equal to 0.9*1/max(eigenvalue)
% @param: Power - 1/0 - If 1, the beta term is flipped to -beta
% @param: G - The matrix on which to calculacte centrality
% @return: centralities - Returns a vector of centralities
%%%
function centralities = BonacichCentrality(alpha,beta,normAlpha,normBeta,Power,G)
n=length(G);

if normBeta==1
    if abs(max(eig(G))) >0
    beta=0.9*abs(1/max(eig(G)));
else % Nilpotent matrix
    beta=1;
    end
    
end
if Power==1
   beta=-beta; 
end

B=1.*((eye(n,n)-beta.*G)^(-1))*G*ones(n,1);

if normAlpha==1
    
    alpha=sqrt(n/(B'*B));
end

centralities = B.*alpha;
end
