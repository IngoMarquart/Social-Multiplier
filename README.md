# What raises the social multiplier? Embedding versus decoupling as formal design choices shaping informal networks


Github repository for the paper

## Abstract

What design-choice should a newly-appointed CEO choose to shape informal networks among her employees, so that their performance is better than if they worked in isolation? What formal intervention, in other words, should she pursue to raise the social multiplier—that part of her organization’s productivity brought forward by social influence? To address this question, we propose an agent-based network model that contrasts two opposing design-choices: embedding and decoupling. Embedding is a formal intervention geared towards raising productivity by increasing peer monitoring and thus ratcheting up social pressure. Conversely, decoupling seeks to curtail social influence and thus reduce distraction. Our analyses reveal four main patterns. First, decoupling is preferable to embedding for an organization whose skill-level distribution is left-skewed and whose employees are more likely to self-enhance (by fixating on less-skilled colleagues) than to self-improve (by eying those with more skill). Second, attempts to rewire informal networks are increasingly irrelevant under what Blau (1977) termed positive consolidation: when skill is concentrated more among self-improvers than among self-enhancers. Third, the mixture of embedding and negative consolidation can catalyze the emergence of two anomalous types of opinion leaders: corrupted role-models—unexpectedly poised to destroy an otherwise well-designed organization; and inspiring underdogs—surprisingly capable of revitalizing an organization otherwise derailed by obsessive status seeking. Fourth, social multipliers are rapidly self-reinforcing, imposing pressure on social planners to act quickly. Our results cast new light on the emergence of attention networks and carry implications for research at the intersection of formal and informal structure.

# Readme


This simulation calculates the subgame perfect Nash equilibrum of the game studied in the paper, under full and bounded rationality. 

Consider options in config.m to set up the simulation. Run it in main.m

Note the following points

- For graphing, parfor loops in main.m should be disabled

- For non-convex utility (v<1), the Matlab global optimization toolpackage is required, and will be used automatically

- Non-convex utility is difficult to compute, so smaller samples are recommended. To reproduce the state space, set v=1.

- Samples used in the paper are given as configuration settings. Please comment and uncomment accordingly in config.m, then run main.m
