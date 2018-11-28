# Comments #

* Good use of Git, commit messages became reasonable over time
* But do make sure you have the same settings for your name on all machines — glancing over stuff in bigger projects would make it confusing to sometimes see Ingo and sometimes IngoMarquart.
* Not sure what happend at 781d31633 vs 4a41b76 -- looks like pretty much the same commit done twice?
* I would not add dates on branches. You see the dates on commits anyhow. Rather focus on the content than on meta-information
* Even some tests, even if failed :-) Make sure to use some automated testing framework, this would be a good place to start https://www.mathworks.com/help/matlab/matlab-unit-test-framework.html
* The code is complex enough to merit separation into various directories. Think the three blocks input (simulation of raw data) -- analysis -- output (tables, graphs)
* 1.3


# SocArch

"Shaping the social architecture of a start-up: What raises the social multiplier?" 
Author: Ingo Marquarts, Nghi Truong, Matthew Bothner, Richard Haynes

This simulation calculates the subgame perfect Nash equilibrum of the game studied in the paper.
mainSPNE.m is the start script, which sets up to generate a large dataset of such equilibria.
It alls SimulateAttSubgame.m with different random seeds and parameters.
