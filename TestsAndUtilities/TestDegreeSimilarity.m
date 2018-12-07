%% Work in Progress
classdef TestDegreeSimilarity < matlab.unittest.TestCase

    properties (TestParameter)
    [0,0,0,0; ...
     0,0,0,0; ...
     0,0,0,0; ...
     0,0,0,0];

    	matrices={'Empty', ...
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
    		'SymFollower', ... 
    		[0,1,0,0; ...
    		 1,0,1,1; ...
    		 0,1,0,0; ...
    		 0,1,0,0], ...
    	}
        type = {'single','double','uint16'};
        level = struct('small', 2,'medium', 4, 'large', 6);
        side = struct('small', 9, 'medium', 81,'large', 729);
    end