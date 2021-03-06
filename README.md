# GlobalOptimality

This page contains codes for numerical solution of globally optimal feedback control protocol for qubit purification in presence of measurement inefficiency.  The codes are all written in matlab.

The codes in this repository were used to generate the solutions in the paper "Optimality of feedback control for qubit purification under inefficient measurement", by Y. Jiang, X. Wang, L. Martin, and K. B. Whaley, arXiv:1910.02487.
Usage of these codes should reference that paper.

There are 4 codes here:

Optimal.m computes the globally optimal control look up table, specifying the value of u for each (r,t). 

run_feedback.m simulates the trajectories when feedback is implemented at all values (r,t).

run_nofeedback.m simulates the trajectories when no feedback is implemented at all values (r,t).

run_Optimal.m simulates the trajectories when the value of u specified by the globally optimal control protocol is used at each value of (r,t). 

