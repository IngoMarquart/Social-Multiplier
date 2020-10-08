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

% New utility
n=length(theta);
if e<1
    ebar=e./(1-e);
    D=eye(n,n).*sum(A,2);
    O=eye(n,n)+ebar.*D;
    F=theta+ebar*A*x_t_1;
    x=O\F;
else
    x=A*x_t_1;
end
% Old utility
% n=length(theta);
% D=eye(n,n).*sum(A,2);
% O=eye(n,n)+e.*D;
% x=O\(theta+e.*A*x_t_1);
%x=(eye(n,n)+e.*(eye(n,n).*sum(A,2)))\(theta+e.*A*x_t_1);


end