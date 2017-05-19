Short term:
- do commit

- run on sherlock

- KKT system S_diag and X_diag as variables

- *first stable step find optimal delta*
- fix delta recording in output
- delta, flags, was the last iteration a failure?
- increase delta when ever there is any sort of failure

- *fix error with predicted reduction of complementarity*

- add option to move dual and primal iterates independently

- mu choice?

- add protections to aggressive steps
- move dual slowly? maybe not i think one is better off regularizing.
#- corrections, only do if predicted progress is good (and stop line search immediately if it is not).
- solve MUMPS issues, version
#- add symmetric KKT system solver and deal with inaccuracy
#- write proper line search for stabilization step
#- non-linear updates of primal variables
- infeasiblity detection
- unboundedness detection
- momentum/homogenous style scaling

- create notes of what I am doing
- validate LP direction

Long term:
- re-use permutations for cholesky
- filter during stable steps, either accept improvement in dual or primal

Blue sky
- run non-linear CG to 0.99 accuracy in 100 iterations. Use NC certificates to deduce correct delta.
- do linear algebra to deal with duplicates of constraints
