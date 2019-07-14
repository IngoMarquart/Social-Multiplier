%%
%   ceoAction.m
%   runs one T on the given firm
%%
% @param: firm- a firm
% @return: firm - a firm
%%
function firm = ceoAction(firm)

% Act whenever half time point is reached
if firm.T>=floor(firm.maxT/2)
    switch firm.ceoAct
        case "Double"
            firm.e=100000;

        case "Half"
            firm.e=firm.e/2;

        case "Zero"
            firm.e=0;

    end
    % We allow for only one CEO action
    firm.startCeoAct=firm.ceoAct;
    firm.ceoAct="Off";
end


end