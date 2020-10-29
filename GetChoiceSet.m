function [ChoiceCell, NrChoices]=GetChoiceSet(G,n)
% GetChoiceSet - For each actor, calculates a Cell containing the indecies of peers who are neighbors in G
%%
% @param: G,n - G network of n actors
% @return: ChoiceCell - n-dimensional cell containing vectors detailing the indecies of neighbors of i
% @return: NrChoices - number of choices for actor i
%%
ChoiceCell={};
NrChoices=zeros(1,n);
for i =1:n
    row=G(i,:);
    ChoiceCell{i}=find(row);
    NrChoices(i)=length(find(row));
end

end