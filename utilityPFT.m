%%
% % utilityPFT
% Returns utility
% @param: x_i - output if i in t
% @param: x_t_1 - last period output
% @param: A - Attention matrix
% @param: e - Embedding
% @param: theta - theta levels
% @param: PsiVec - Vector of Psi benefits
% @param: maxDegree - Maximum number of peers
% @param: conParam - Convex or concave influence of monitoring
% @return: u - Negative Utility Value
%%

function u=utilityPFT(x_i,x_t_1,a_i,theta,e,PsiVec,i,maxDegree,conParam)

    x_t_1=x_t_1(:);
    PsiVec=PsiVec(:);
    x=x_t_1;
    x(i)=x_i;
    a_i=a_i(:);
    theta=theta(:);
    n=length(theta);
    ez=ones(n,1);
    % Private utility, positive part
    PrivUtil=(x(i)-theta(i))^2;
    % Make sure no connection to oneself
    PsiVec(i)=-100;
    % Expected benefit
    if sum(a_i>0) > maxDegree
        CBenUtil=-100;
    elseif conParam==0
        CBenUtil=((a_i'.^(1))*(PsiVec));
    else
        CBenUtil=((a_i'.^(conParam))*(PsiVec));
    end
    
    % Expected non-alignment cost
    ConfUtil=a_i'*((x(i).*ez-x).*(x(i).*ez-x));
    if sum(a_i)>0
    % Full utility new model
        u=-PrivUtil.*(1-e)+e.*CBenUtil-e.*ConfUtil;
    else
        u=-0;
    end
end