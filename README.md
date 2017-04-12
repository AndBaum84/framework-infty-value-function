# framework-infty-value-function
In this branch a numerical framework is established which can be applied to investigate the dynamic, capital-theoretic model evolved by Board and ter Vehn (2013). Solving the problem of the computation of the infinite-horizon value functions, coupled Hamilton-Jacobi-Bellman equations, are treated to build systems of linear equations under optimal control. The resulting value-curves are approximated by Chebyshev-polynomials to receive a useful form of the value functions for further analysis. In combination with the definition of an equilibrium in the model, Newton's method in combination with the bisection-method are adopted to find optimal values. This procedure is applied to different system-settings to evaluate potential boundaries of the model.

The results of this approach provide evidence that the marginal flow cost of an investment decision and the discount rate, restrain the existence of an equilibrium in the same direction. More generally, the existence of an equilibrium is depending on the assessed values for the system-parameters and their interrelation. Depending on the stationary point of the reputational dynamics, more than one equilibrium exists for certain system-settings. The results are investigated according to their robustness, to proof the validity of the numerical approach. The numerical framework applied in this paper is a useful approach to identify the unknown value functions of the model Reputation for Quality which can serve as a basis for further research.

1) Routine.m the main schript is: 
	- That scribt will execute the Algorithm
	- All generated data, hence candidates and equilibrium cutoff-values are written into matrices
	- Furthermore, the corresponding Chebyshev-approximation errors, information about the monotonicity and system-parameters are saved as well during the routine

2) plotvaluecurves.m:
	- Plots the value curves according to a cutoff-value and system-parameters that has to be prespecified.

3) printcandidates.m, printoptimumcflex.m, printoptimumrflex.m, printvalueofqualigrid.m:
	- Plot different graphics with output from Routine.m
