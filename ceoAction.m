%%
%   ceoAction.m
%   runs one T on the given firm
%%
% @param: firm- a firm
% @return: firm - a firm
%%
function firm = ceoAction(firm)
if firm.T>2
    diffM=firm.diffM(firm.T-1);
    firm.eMat(firm.T)=max(0.1,firm.eMat(firm.T-1)+sign(diffM));
end


end