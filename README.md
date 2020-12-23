# Codes of "Global dynamics in nonconservative nonlinear Schroedinger equations"

This repository contains the MATLAB codes associated with the paper:
"Global dynamics in nonconservative nonlinear Schr\"odinger equations"
by J Jaquette, J-P Lessard and A Takayasu. ([arXiv:2012.09734 [math.AP]](https://arxiv.org/abs/2012.09734))

**Abstract**  In this paper, we study the global dynamics of a class of nonlinear Schr\"odinger equations using perturbative and non-perturbative methods. We prove the semi-global existence of solutions for initial conditions close to constant. That is, solutions will exist for all positive time or all negative time.  The existence of an open set of initial data which limits to zero in both forward and backward time is also demonstrated. This result in turn forces the non-existence of any real-analytic conserved quantities. For the quadratic case, we prove the existence of two (infinite) families of nontrivial unstable equilibria and prove the existence of heteroclinic orbits limiting to the nontrivial equilibria in backward time and to zero in forward time. By a time reversal argument, we also obtain heteroclinic orbits limiting to the nontrivial equilibria in forward time and to zero in backward time. The proofs for the quadratic equation are computer-assisted and rely on three separate ingredients: an enclosure of a local unstable manifold at the equilibria, a rigorous integration of the flow (starting from the unstable manifold) and a proof that the solution enters a validated stable set (hence showing convergence to zero).

These codes require *MATLAB* with [*INTLAB* - INTerval LABoratory](http://www.ti3.tu-harburg.de/rump/intlab/) (MATLAB toolbox for interval arithmetic) version 11 and [*Chebfun* - numerical computing with functions](https://www.chebfun.org/) version 5.7.0.

---

A rough correspondence for some of the files & computational procedures in the paper are as follows:

### Existence of a steady state and an eigenpar for NLS

```
>> cd Proofs_Eigenpairs
>> script_verify_eigenpairs
```

### Constructing a part of unstable manifold via Parameterization method

```
>> cd ../Manifolds/
>> script_get_manifold
```

### Rigorous integration of a flow and validating semi-global existence limiting to zero

```
>> cd ../verify_solution/
>> script_proof_NLS_from_P_at_1
```

Then it proves the existence of heteroclinic orbits limiting to the nontrivial equilibria in backward time and to zero in forward time.


Copyright (C) 2020  J Jaquette, J-P Lessard and A Takayasu.
