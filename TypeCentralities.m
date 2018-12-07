%%%
% % TypeCentralities 
% Return a couple of measures that provide risk
% @param: G - The matrix on which to calculacte centrality
% @param: identities - vector of types in [-1,1]^n
% @return: centralities - Returns a vector of centralities
%%%
function centralities = TypeCentralities(G, identities)
n=length(G);
if G==G'
	Gsym=True
end

%% 1. Out Centrality
outcentrality=BonacichCentrality(0,0,1,1,0,G)
outcentrality=outcentrality.*identities;

%% 2. In Centrality
incentrality=BonacichCentrality(0,0,1,1,0,G)
incentrality=incentrality.*identities;

%% 3. Influence Centrality
influenceG=repmat(identities,n,1)*G;
influencecentrality=1.*((eye(n,n)-beta.*G)^(-1))*influenceG*ones(n,1);
alpha=sqrt(n/(influencecentrality'*influencecentrality));
influencecentrality = influencecentrality.*alpha;


centralities.outcentrality=outcentrality
centralities.incentrality=incentrality
centralities.influencecentrality=influencecentrality
end