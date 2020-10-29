function a_i = RecoverPi(aiRest, ChoiceCell, n)
% RecoverPi - From a restricted a vector (in the space of valid choices), recover a A-matrix conform vector (in the space of n actors)
%%
% @param: n - number of actors
% @param: ChoiceCell - Vector containing the indecies of neighbors
% @param: aiRest - Chosen weights among neighbors
% @return: a_i - Chosen weights among neighbors and non-neighbors in G
%% 
    a_i = zeros(n,1);
    a_i(ChoiceCell)=aiRest;
end