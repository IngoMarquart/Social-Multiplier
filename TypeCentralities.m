%%%
% % TypeCentralities
% Return a couple of measures that provide risk
% @param: G - The matrix on which to calculacte centrality
% @param: identities - vector of types in [-1,1]^n
% @return: centralities - Returns a vector of centralities
%%%
function centralities = TypeCentralities(G, identities)
n=length(identities);
if rank(G)==0
    centralities.outcentrality=zeros(n,1);
    centralities.incentrality=zeros(n,1);
    centralities.influencecentrality=zeros(n,1);
    %warning("G has zero rank")
else
    if abs(max(eig(G))) >0
        beta=0.9*abs(1/max(eig(G)));
    else % Nilpotent matrix
        beta=1;
    end
    %% 1. Out Centrality
    outcentrality=BonacichCentrality(0,0,1,1,0,G);
    outcentrality=outcentrality.*identities;
    
    %% 2. In Centrality
    incentrality=BonacichCentrality(0,0,1,1,0,G');
    incentrality=incentrality.*identities;
    
    %% 3. Influence Centrality
    influenceG=repmat(identities,1,n).*G';
    influencecentrality=1.*((eye(n,n)-beta.*G')^(-1))*influenceG*ones(n,1);
    if influencecentrality'*influencecentrality == 0
        alpha=0;
    else
        alpha=sqrt(n/(influencecentrality'*influencecentrality));
    end
    influencecentrality = influencecentrality.*alpha;
    
    centralities.outcentrality=outcentrality;
    centralities.incentrality=incentrality;
    centralities.influencecentrality=influencecentrality;
end
end