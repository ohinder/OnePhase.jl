Short term:
- cutest
- increase delta when ever there is any sort of failure
- move dual slowly?
#- corrections, only do if predicted progress is good (and stop line search immediately if it is not).
- solve MUMPS issues, version
#- add symmetric KKT system solver and deal with inaccuracy
#- write proper line search for stabilization step
- non-linear updates of primal variables
- infeasiblity detection
- unboundedness detection
- momentum/homogenous style scaling

- create notes of what I am doing
- validate LP direction

Long term:
- re-use permutations for cholesky
- filter during stable steps, either accept improvment in dual or primal

Blue sky
- run non-linear CG to 0.99 accuracy in 100 iterations. Use NC certificates to deduce correct delta.
- do linear algebra to deal with duplicates of constraints
