# SocArch

"Shaping the social architecture of a start-up: What raises the social multiplier?" 
Author: Ingo Marquarts, Nghi Truong, Matthew Bothner, Richard Haynes

This simulation calculates the subgame perfect Nash equilibrum of the game studied in the paper.
mainSPNE.m is the start script, which sets up to generate a large dataset of such equilibria.
It alls SimulateAttSubgame.m with different random seeds and parameters.

The variable names do not reflect current changes in the paper.

It was attempted to refactor all functionality into its own matlab .m files. This means the natural way to read this code is by starting in mainSPNE.m, for the grouped runs, and then SimulateAttSubgame.m, for each single run of a firm.

# Optimization algorithms

Note how there are different ways to derive the equilibrium, global, local and discrete, in ascending order of speed. Global requires the matlab global optimization toolbox. It is robust, and finds the global optimum, but not fast enough to calculate the whole state space. 

Local is faster, but we can not formally guarantee that the equilibrium wrt. to attention is unique. Tests indicate that the local optimization algorithm fails if n<10, thus a hard-coded fallback is included.

Both algorithms are run purely numerically without further inputs such as the Hessian. This is on purpose, as we want the simulation to validate our analytical results from first-order assumptions.

The discrete algorithm is very fast. It is based on results we derive in the model, in particular that peer choice is unique, leading to a simple linear algebra & decision problem. 