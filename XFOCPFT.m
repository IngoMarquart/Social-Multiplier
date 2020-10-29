function x=XFOCPFT(x_t_1,A,theta,e,k)
% XFOCPFT - Solves FOC for all actors, returning PFT output
%%
% @param: x_t_1 - Prior period output
% @param: A - Attention matrix
% @param: e - Embedding
% @param: theta - theta levels
% @return: x - vector of effort levels
%%
% New utility
n=length(theta);

x_bar=A*x_t_1;

logic_vec=theta<=x_bar;

% Left leg of utility function
if e<1
    ebar=e./(1-e);
    D=eye(n,n).*sum(A,2);
    O=eye(n,n)+ebar.*D;
    F=theta+ebar*A*x_t_1;
    x_n=O\F;
else
    x_n=A*x_t_1;
end

% Right leg of utility function
if e<1
    ebar=(e*k)./(1-e);
    D=eye(n,n).*sum(A,2);
    O=eye(n,n)+ebar.*D;
    F=theta+ebar*A*x_t_1;
    x_k=O\F;
else
    x_k=A*x_t_1;
end

% Assing left or right leg of utility
x=logic_vec.*x_n+(1-logic_vec).*x_k;


end
