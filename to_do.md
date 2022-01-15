- low tolerance Ipopt
- infeasible Ipopt

Short term:
- mu autotune scaling -- run on whole test set.
- increase delta to min stable line search??? *test*
- new delta strategy with loglog rate???
- write code that saves all the CUTEst files need for the paper in one go! (i.e., have an experimental script and a paper producing script)

- add COPS problems
- improved initial delta selection

Medium term:
- detect is problem is caused by factorization issues or lack of smoothness of functions (identify function, a direction and a point)
- *improve efficiency of schur complement and eval_jac*

- deal correctly with exceptions
- add parameters for termination criterion
- test unconstrained problems

*- create infeasible test set*
*- scaled termination criterion*
*- run full netlib test*

*- mu stuff *

*- initialization scheme*
- *fix error with predicted reduction of complementarity*
- **stabilization should prioritize complementarity if dual is small and comp not satisfied**
- *measure complementarity in output in relative terms*

- split up parameters i.e. initialization parameters etc ...

- install IPOPT on sherlock


- KKT system S_diag and X_diag as variables

- increase delta when ever there is any sort of failure


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

**Long term:**
- (1) find LP solution first, (2) start from analytic centre, (3) re-write so problem is well-conditioned

- automatic scaling trust region algorithm

- momentum/CG in stabilization steps

- re-use permutations for cholesky
- filter during stable steps, either accept improvement in dual or primal

**Blue sky**
- run non-linear CG to 0.99 accuracy in 100 iterations. Use NC certificates to deduce correct delta.
- do linear algebra to deal with duplicates of constraints
