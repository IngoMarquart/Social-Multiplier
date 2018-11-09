%%
% % JacksonRogersNW
% Creates a G Network following Jackson&Rogers(2014) algorithm
% @param: n - Number of actors
% @param: mn - Number of peers to find in either of the two steps
% @param: pn - Probability to connect to a found peer
% @param: randomseed - Randomseed to set.
%% 

function G=JacksonRogersNW(n,mn,pn, randomseed)

    % Set random seed 
    s = RandStream('mcg16807','Seed',randomseed);
    RandStream.setGlobalStream(s);
    % Assumption in Jackson and Rogers, pr=pn
    pr=pn;
    % Set n-based default values
    if mn==0
        mr=min(round(n/12),3);
        mn=min(round(n/12),3);
    else
        mr=mn;
    end


    % We will add starting nodes that we will delete later
    nfix=n+mr+mn+1;
    t=mr+mn+1;

    % Create fully connected network
    G=ones(t,t)-eye(t,t);

    for period = t:nfix-2
        oldG=G;
        % Step 1: select uniformly at random
        selectvec=randsample(t,mr);
        row=zeros(1,t+1);
        % With probability pr, connect to those nodes
        for link = selectvec(:)'
            row(link) = 1.*binornd(1,pr);
        end
        % Step 2: Collect a set of nodes connected to parents
        parents = selectvec;
        set=[];
        % For each parent node
        for i = parents(:)'
            % Collect links of parents
            prow=oldG(i,:);
            % Collect indeces of those nodes and merge with parent set
            set = union(find(prow),set);
        end
        % Remove parents from this set
        set = setdiff(set,parents);
        
        rown=zeros(1,t+1);
        % With probability pn, connect to those nodes
        for link = set(:)'
            rown(link) = 1.*binornd(1,pn);
        end
        
        % Final connections: add rows
        row=row+rown;
        
        % Construct new matrix and proceed one period
        t=t+1;
        G=zeros(t,t);
        G(1:t-1,1:t-1)=oldG;
        G(t,:)=row;
        G(:,t)=row';
    end
    G=G(1:n,1:n);
    % % Remove starting nodes
    % G=G(mr+mn+1:end,mr+mn+1:end);
    %
    % GPgraph=graph(G);
    % h=plot(GPgraph,'Layout','force')
end