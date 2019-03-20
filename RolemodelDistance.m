%%%
% % RolemodelDistance
% This function returns the shortest path (in a symmetric network sense) between most central actors
% @param: G - The matrix returning distances
% @return: centralities - Returns a vector of centralities
%%%
function calcdistance = RolemodelDistance(G)
	
	% Symmetrize Network
	G=(G+G')./2;
	G(G>0)=1;
    G(G<0)=0;
	% Sort centralities
	Bcent=BonacichCentrality(0,0,1,1,0,G);
	% Sort the centralities, G and create Graph
	[Bcent,sortIndex]=sort(Bcent,'descend');
	GGraph=graph(G);
    
	[~,d] = shortestpath(GGraph,sortIndex(1),sortIndex(2));
	calcdistance=max(d,0);
end