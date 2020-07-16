%%
% DiscreteChoice
% Assuming that each actor picks only one peer, this function find the one
% giving the most utility
% @param: a_t_1 - Attention matrix of last time period
% @param: e,delta - Embedding and scaling variable
% @param: theta_i, theta - theta vector and theta of i
% @param: Psi - Expected benefit function handle
% @param: Choice - Vector given the indecies of valid peers (where connection in G exists)
% @param: nrChoice - Number of choices
% @return: util - Utility value
% @return: a_i_star - Attention choice in the space of all n actors
%%
function [util,a_i_star]=ConvexDiscreteChoicePFT(a_t_1,x_t_1,e, theta, Psi,i,Choice,nrChoice, rationality, conParam)

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
    % FOC Version: x anticipated based on a_t_1
    if rationality > 0
        x_rat=XFOCSPNE(a,theta,e);
        x=(1-rationality).*x_pft+rationality.*x_rat;
    else
        x=x_pft;
    end
    
    % Private utility, positive part
    PrivUtil=(x(i)-theta(i))^2;
    % Calculate a benefit vector for each potential peer
    PsiVec = Psi(theta(i),theta);
    % Make sure no connection to oneself
    PsiVec(i)=-100;
    % Expected benefit
    if sum(a_i>0) > 10
        CBenUtil=0;
    elseif conParam==0
        CBenUtil=((a_i'.^(1))*(PsiVec));
    else
        CBenUtil=((a_i'.^(conParam))*(PsiVec));
    end
    
    % Expected non-alignment cost
    ConfUtil=a_i'*((x(i).*ez-x).*(x(i).*ez-x));
    
    % Full utility new model
    NewUtil=-PrivUtil.*(1-e)+e.*CBenUtil-e.*ConfUtil;
    
    % Full utility old model
    %NewUtil=-PrivUtil+e.*CBenUtil-e.*ConfUtil;
    
    % If this choice is better than the previous one, do this
    if NewUtil > util
        a_i_star=a_i;
        util=NewUtil;
    end
end






end