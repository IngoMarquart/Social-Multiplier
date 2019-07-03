%%
% % TaskNetwork
% Creates a G Network following Jackson&Rogers(2014) algorithm
% @param: n - Number of actors
% @param: mn - Number of peers to find in either of the two steps
% @param: pn - Probability to connect to a found peer
% @param: randomseed - Randomseed to set.
%% 

function [G,cluster,links]=TaskNetwork(n,cluster,modularity, symmetry, randomseed)
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
% Collect zero
zeroG=find(G==0);
numZeroG=numel(find(G==0));

if modularity >0
    % modularity -> percentage of links that should be nonzero
    links=floor(numZeroG*modularity);
    if symmetry==0
        links=floor(links./2);
    end
    % Index those onto zero elements
    gIndex=zeroG(randperm(numZeroG,links));
    % Set G to 1
    G(gIndex)=1;
    if symmetry==1
        G=G+G';
        G=G>0;
    end
end

% Delete diagonals
G=G-eye(n,n);
G=G==1; % May create negative links for non clustered nodes
end