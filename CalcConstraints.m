function [ConA,Conb] = CalcConstraints(n, NrChoices)
ConA={};
Conb={};
%% Constraints
for i = 1:n        
    
        % Build a fake prest_prev vector based on choice sets
        prest_prev = zeros(NrChoices(i),1);
        % Constraint matrices determine
        % 1. Each variable less or equal 1
        % 2. Each variable greater or equal 0
        % 3. The sum is smaller or equal 1
        % 4. The sum is larger or equal 0
        ConA{i}=[eye(length(prest_prev));-1.*eye(length(prest_prev));ones(1,length(prest_prev));-1.*ones(1,length(prest_prev))];
        Conb{i}=[ones(length(prest_prev),1);zeros(length(prest_prev),1);1;0];
end       
        
end