%%
% % TaskNetworkAss
% Creates a G Network representing a cluster-wise assembly line
% @param: n - Number of actors
% @param: cluster - number of clusters
% @param: modularity - percentage of off-diagonal links to be 1
% @param: fwRate - Percentage of links forward connection between clusters
% @param: randomseed - Randomseed to set.
%%

function [G,cluster,links]=TaskNetworkAss(n,cluster,modularity, assembly, gSymmetry, randomseed)

if gSymmetry==1
    fwRate=0.5;
else 
    fwRate=1;
end

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

% Set non clustered diagonals to 1
% Delete diagonals
G=G+eye(n,n);
G=G>=1;

% How many links in the diagonal
numOneG=numel(find(G==1));

if modularity >0
    % modularity -> percentage of links that should be nonzero
    links=floor(numOneG*modularity);
    linksAss=round(links.*assembly);
    linksRan=links-linksAss;
    %% Assembly line process
    % How many links should be to forward vs. backward direction
    linksFW=round(linksAss.*fwRate);
    linksBW=round(linksAss.*(1-fwRate));
    % Get linear indecies
    fillFW=find(G==1);
    % Permute forward links and backward links
    gIndexFW=fillFW(randperm(numOneG,linksFW));
    gIndexBW=fillFW(randperm(numOneG,linksBW));
    % Create two diagonal matrices with fw and bw links set to 1
    fwMat=zeros(n,n);
    bwMat=zeros(n,n);
    fwMat(gIndexFW)=1;
    bwMat(gIndexBW)=1;
    % Circular shift both matrices off the diagonal
    % Ignoring the nodes not in clusters (cRest)
    fwMat(1:end-cRest,:)=circshift(fwMat(1:end-cRest,:),cSize,2);
    bwMat(:,1:end-cRest)=circshift(bwMat(:,1:end-cRest),cSize,1);
    % Delete overlap
    fwMat(n-cSize:n-cRest,1:cSize-cRest)=0;
    bwMat(1:cSize-cRest,n-cSize:n-cRest)=0;
    % Delete diagonal elements of non-clustered nodes
    fwMat=fwMat-eye(n,n).*fwMat;
    bwMat=bwMat-eye(n,n).*bwMat;
    % Add both matrices
    G=G+fwMat+bwMat;
    
    %% Random Links
    % Collect zero
    zeroG=find(G==0);
    numZeroG=numel(find(G==0));
        % Index those onto zero elements
    gIndex=zeroG(randperm(numZeroG,linksRan));
    % Set G to 1
    G(gIndex)=1;
end

% Delete diagonals
G=G-eye(n,n).*G;
G=G==1; % May create negative links for non clustered nodes
end