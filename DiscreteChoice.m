%%
% DiscreteChoice
% Assuming that each actor picks only one peer, this function find the one
% giving the most utility
% @param: P_t_1 - Attention matrix of last time period
% @param: g,delta - Embedding and scaling variable
% @param: theta_i, theta - theta vector and theta of i
% @param: Psi - Expected benefit function handle
% @param: Choice - Vector given the indecies of valid peers (where connection in G exists)
% @param: nrChoice - Number of choices
% @return: util - Utility value
% @return: p_i_star - Attention choice in the space of all n actors
%%
function [util,p_i_star]=DiscreteChoice(P_t_1,g,delta,theta_i, theta, Psi,i,Choice,nrChoice, rationality)

theta=theta(:);
n=length(theta);
% Define One vector
ez=ones(n,1);

% Default values
util=0;
p_i_star=zeros(n,1);

% Iterate over choice sets
for z = 1:nrChoice
    % Create a prest
    prest=zeros(nrChoice,1);
    prest(z)=1;
    % Recover p_i from p restriction
    p_i = RecoverPi(prest, Choice, n);
    ebar=g/(g+1);
    % Create new P matrix
    P=P_t_1;
    % Assemble new P
    P(i,:)=p_i';
    % Get expected x
    x=XFOCSPNE(P,delta,theta,g);
    x_i=x(i);
    % Calculate boundedly rational or rational choice
    x=(1-rationality).*theta+rationality.*x;
    %x(i)=x_i;
    x(i)=(1-ebar).*theta_i+ebar.*p_i'*x;

    % Private utility, positive part
    PrivUtil=(x_i-theta_i)^2;
    % Calculate a benefit vector for each potential peer
    PsiVec = Psi(theta_i,theta);
    % Make sure no connection to oneself
    PsiVec(i)=-100;
    % Expected benefit
    CBenUtil=p_i'*(PsiVec);
    % Expected non-alignment cost
    ConfUtil=p_i'*((x_i.*ez-x).*(x_i.*ez-x));
    % Full utility
    NewUtil=-PrivUtil+g^(1).*delta.*CBenUtil-g.*delta.*ConfUtil;
    % If this choice is better than the previous one, do this
    if NewUtil > util
       p_i_star=p_i;
       util=NewUtil;
    end
end






end