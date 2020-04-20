%%
% % XFOCPFT
% Solves FOC for all actors, returning PFT output
% @param: x_t_1 - Prior period output
% @param: A - Attention matrix
% @param: e - Embedding 
% @param: theta - theta levels
% @return: x - vector of effort levels
%% 

function x=XFOCPFT(x_t_1,A,theta,e)
    theta=theta(:); % Make sure column vector
    x_t_1=x_t_1(:); % Make sure column vector
    x=(1-e)*theta+e*A*x_t_1;
    
end