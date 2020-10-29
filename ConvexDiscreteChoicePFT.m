function [util,a_i_star]=ConvexDiscreteChoicePFT(a_t_1,x_t_1,e, theta, Psi,i,Choice,nrChoice, rationality, conParam, maxDegree)
% ConvexDiscreteChoicePFT - Assuming that each actor picks only one peer, this function find the one giving the most utility
%%
% @param: a_t_1 - Attention matrix of last time period
% @param: e,delta - Embedding and scaling variable
% @param: theta_i, theta - theta vector and theta of i
% @param: Psi - Expected benefit function handle
% @param: Choice - Vector given the indecies of valid peers (where connection in G exists)
% @param: nrChoice - Number of choices
% @return: util - Negative Utility value
% @return: a_i_star - Attention choice in the space of all n actors
%%
theta=theta(:);
n=length(theta);
% Define One vector
ez=ones(n,1);

% Default values
util=0;
a_i_star=zeros(n,1);

% Iterate over choice sets
for z = 1:nrChoice
    
    % Create a prest
    prest=zeros(nrChoice,1);
    prest(z)=1;
    % Recover p_i from p restriction
    a_i = RecoverPi(prest, Choice, n);
    % Create new P matrix
    a=a_t_1;
    % Assemble new P
    a(i,:)=a_i';
    % Get expected x
    
    % PFT Version: x_t_1 for j != i
    x_pft=x_t_1;
    x_pft_foc=XFOCPFT(x_t_1,a,theta,e);
    x_pft(i)=x_pft_foc(i);
    x=x_pft;
    

    % Calculate a benefit vector for each potential peer
    PsiVec = Psi(theta(i),theta);

    % Get utility value
    NewUtil=-utilityPFT(x(i),x,a_i,theta,e,PsiVec,i,maxDegree,conParam);
    
    % If this choice is better than the previous one, do this
    if NewUtil < util
        a_i_star=a_i;
        util=NewUtil;
    end
end






end