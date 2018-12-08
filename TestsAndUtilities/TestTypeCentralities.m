%%%
% TEST
% % TypeCentralities
% Return a couple of measures that provide risk
% @param: G - The matrix on which to calculacte centrality
% @param: identities - vector of types in [-1,1]^n
% @return: centralities - Returns a vector of centralities
%%%
classdef TestTypeCentralities < matlab.unittest.TestCase
    properties (ClassSetupParameter)
        identities = {[1,1,1,1]',[0,0,0,0]',[1,-1,1,-1]',[-1,1,-1,1]',[-1,-1,-1,-1]'};
    end
    properties (TestParameter)
        matrices=struct('Empty', ...
            [0,0,0,0; ...
            0,0,0,0; ...
            0,0,0,0; ...
            0,0,0,0], ...
            'SymPairs', ...
            [0,1,0,0; ...
            1,0,0,0; ...
            0,0,0,1; ...
            0,0,1,0], ...
            'SymCentral', ...
            [0,1,1,1; ...
            1,0,0,0; ...
            1,0,0,0; ...
            1,0,0,0], ...
            'SymFollower', ...
            [0,1,0,0; ...
            1,0,1,1; ...
            0,1,0,0; ...
            0,1,0,0], ...
            'Central', ...
            [0,0,0,0; ...
            1,0,0,0; ...
            1,0,0,0; ...
            1,0,0,0], ...
            'Follower', ...
            [0,1,0,0; ...
            0,0,0,0; ...
            0,1,0,0; ...
            0,1,0,0])
        typecents=struct('Empty', [1,1] , ...
            'SymPairs', [1,1] , ...
            'SymCentral', [1,1], ...
            'SymFollower', [1,1],...
            'Central', [1,0],...
            'Follower', [0,1])
    end
    methods (TestClassSetup)
        function ClassSetup(testCase, identities)
            identity=identities;
        end
    end
    methods (Test)
        function testSignage(testCase, identities, matrices)
            centralities=TypeCentralities(matrices, identities);
            if identities(1)>=0
                testCase.verifyGreaterThanOrEqual(centralities.outcentrality(1),0);
                testCase.verifyGreaterThanOrEqual(centralities.incentrality(1),0);
            else
                testCase.verifyLessThanOrEqual(centralities.outcentrality(1),0);
                testCase.verifyLessThanOrEqual(centralities.incentrality(1),0);
            end
        end
    end
    methods (Test, ParameterCombination='sequential')
        function testInfluenceCentDirection(testCase, matrices, identities, typecents)
            centralities=TypeCentralities(matrices, identities);
            typecents(1)=typecents(1).*identities(1);
            typecents(2)=typecents(2).*identities(2);
            if typecents(1) ~= typecents(2) && rank(matrices) > 0
                assertion=typecents(1)>typecents(2);
                proposition=centralities.influencecentrality(1)>=centralities.influencecentrality(2);
                
                message=strcat("Assertion: ",num2str(typecents(1))," > ", num2str(typecents(2)));
                message=strcat(message, " equals",num2str(centralities.influencecentrality(1))," > ", num2str(centralities.influencecentrality(2)));
                
                testCase.verifyEqual(assertion, proposition,message);
            end
        end
    end
end
