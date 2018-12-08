%%%
% % DegreeSimilarity 
% Returns a matrix of similarities in degree
% @param: G - The matrix on which to calculacte centrality
% @param: types - vector of types
% @return: centralities - Returns a vector of centralities
%%%
function [outdegreesim,indegreesim] = DegreeSimilarity(G)
n=length(G);


inDeg=G'*ones(n,1);
outDeg=G*ones(n,1);

indegreesim=repmat(inDeg,1,n)-repmat(inDeg',n,1);
outdegreesim=repmat(outDeg,1,n)-repmat(outDeg',n,1);

end