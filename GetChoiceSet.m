function [ChoiceCell, NrChoices]=GetChoiceSet(G,n)

ChoiceCell={};
NrChoices=zeros(1,n);
for i =1:n
    row=G(i,:);
    ChoiceCell{i}=find(row);
    NrChoices(i)=length(find(row));
end

end