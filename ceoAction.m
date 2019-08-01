%%
%   ceoAction.m
%   runs one T on the given firm
%%
% @param: firm- a firm
% @return: firm - a firm
%%
function firm = ceoAction(firm)

% Act whenever half time point is reached
if firm.T>=firm.ceoStartT
    switch firm.ceoAct
        case "Embed"
            firm.e=100000;
            firm.startCeoAct=firm.ceoAct;
        case "LowEmbed"
            firm.e=2;
            firm.startCeoAct=firm.ceoAct;
        case "Decouple"
            firm.e=0;
            firm.startCeoAct=firm.ceoAct;
    end
    % We allow for only one CEO action
    
    firm.ceoAct="Off";
end


end