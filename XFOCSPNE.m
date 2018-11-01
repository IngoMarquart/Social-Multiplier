function x=XFOCSPNE(P,delta,theta,g)

    % Create Laplacian matrix
    n=length(theta);
    DP=eye(n,n).*sum(P,2);
    L=DP-P;
    
    % Calculate X-vector given P
    %x=(eye(n,n)+delta.*L)\theta;
    % Use this BR function if g is in front of conf-cost
    x=(eye(n,n)+g.*delta.*L)\theta;
    
end