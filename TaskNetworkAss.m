%%
% % TaskNetworkAss
% Creates a G Network representing a cluster-wise assembly line
% @param: n - Number of actors
% @param: cluster - number of clusters
% @param: modularity - percentage of off-diagonal links to be 1
% @param: fwRate - Percentage of links forward connection between clusters
% @param: randomseed - Randomseed to set.
%% 

function [G,cluster,links]=TaskNetworkAss(n,cluster,modularity, fwRate, randomseed)
  % Set random seed 
    s = RandStream('mcg16807','Seed',randomseed);
    RandStream.setGlobalStream(s);
    
% We need as least as many agents as cluster
if n<cluster
    cluster=n;
end
% Calculate cluster size
cSize=floor(n./cluster);
cRest=mod(n,cluster);

% Create a network with densely connected clusters
diagG=kron(eye(cluster),ones(cSize));
G=zeros(n);
G(1:(n-cRest),1:(n-cRest))=diagG;
% How many links in the diagonal
numOneG=numel(find(G==1));

if modularity >0
    % modularity -> percentage of links that should be nonzero
    links=floor(numOneG*modularity);
    % How many links should be to forward vs. backward direction
    linksFW=round(links.*fwRate);
    linksBW=round(links.*(1-fwRate));
    % Define diagonal matrix to fill with links
    diagG=kron(eye(cluster),ones(cSize));
    % Get linear indecies
    fillFW=find(diagG==1)
    % Permute forward links and backward links
    gIndexFW=fillFW(randperm(numZeroG,linksFW));
    gIndexBW=fillFW(randperm(numZeroG,linksBW));
    % Create two diagonal matrices with fw and bw links set to 1
    fwMat=zeros(n,n);
    bwMat=zeros(n,n);
    fwMat(gIndexFW)=1;
    bwMat(gIndexBW)=1;
    % Circular shift both matrices off the diagonal
    fwMat=circshift(fwMat,cSize,1);
    bwMat=circshift(bwMat,cSize,1);
    % Delete overlap
    fwMat(n-cSize:n,1:cSize)=0;
    bwMat(1:cSize,n-cSize)=0;
    % Add both matrices
    G=G+fwMat+bwMat;
end

% Delete diagonals
G=G-eye(n,n);
G=G==1; % May create negative links for non clustered nodes
end