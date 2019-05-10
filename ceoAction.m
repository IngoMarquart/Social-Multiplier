%%
%   ceoAction.m
%   runs one T on the given firm
%%
% @param: firm- a firm
% @return: firm - a firm
%%
function firm = ceoAction(firm)
if firm.T>2
    window=1;
    dr=1;
    diffM=firm.diffM(firm.T-1);
    ebar=firm.eMat(firm.T-1)./(firm.eMat(firm.T-1)+1);
    
    
    diffM=firm.diffM;
    if dr==1
    diffM(1)=diffM(2);
     if firm.T<=window
         diffMIndicator=mean(diffM(1:firm.T-1));
     else
         diffMIndicator=mean(diffM((firm.T-window):(firm.T-1)));
     end
    else
    diffMIndicator=diffM(firm.T-1)-diffM(firm.T-2);
    end
    
    new_ebar=min(0.99,max(0,(ebar)+(min(ebar,(1-ebar))./2)*sign(diffMIndicator)));
    firm.eMat(firm.T)=-new_ebar./(new_ebar-1);
    %firm.eMat(firm.T)=max(0.01,(firm.eMat(firm.T-1)+((1-ebar)./1.01)*sign(diffMIndicator)));
end


end