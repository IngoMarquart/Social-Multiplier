%%
% % RecoverPi
% From a restricted p vector (in the space of valid choices), recover a P-matrix conform 
% p vector (in the space of n actors)
% @param: n - number of actors
% @param: ChoiceCell - Vector containing the indecies of neighbors
% @param: prest - Chosen weights among neighbors
% @return: p_i - Chosen weights among neighbors and non-neighbors in G
%% 

function p_i = RecoverPi(prest, ChoiceCell, n)
    p_i = zeros(n,1);
    p_i(ChoiceCell)=prest;
    
    

end