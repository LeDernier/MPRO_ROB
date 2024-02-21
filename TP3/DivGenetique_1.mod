param N;	# first index, nombre d'individus
param Nm; 			# nombre d'individus mâles
param Nf;				# nombre d'individus femelles
param P;	# second index, nombre de paire de chromosomes
param G;	# third index, locus (1e coordonnée du gène dans le chromosome)
param C;	# fourth index, 1er ou 2e chromosome (2e coordonnée du gène dans le chromosome)
param T;				# number of thetas (linear approx of log(.))
param init default 0.001;				# theta_1
param x_UB default 2;	# upper bound on the number of childs


		
param individu{1..N, 1..P, 1..G, 1..C};

## useful parameters and sets generated
set A := {1,2}; # ensemble d'allèles: 1 pour l'allèle dominant ('A') et 2 pour l'allèle recessif ('a')

param proba{n in 1..N, p in 1..P, g in 1..G, a in A} default 1 - 
			((if individu[n,p,g,1] == a then 1 else 0) + 
			(if individu[n,p,g,2] == a then 1 else 0))/2; # probability of non-transmition of the gen
param theta{r in 1..T} := init**((T-r)/(T-1));			  # parameters in the approximation

set I_One{p in 1..P, g in 1..G, a in A} := {n in 1..N: proba[n,p,g,a] = 1};
set I_OneHalf{p in 1..P, g in 1..G, a in A} := {n in 1..N: proba[n,p,g,a] = 0.5};

## decision variables
var x{n in 1..N} integer, >= 0;						# number of childs of person i
var y{p in 1..P, g in 1..G, a in A} >= 0, <= 1; # probability of disappearance of the gen a(p,g,a)
var z{p in 1..P, g in 1..G, a in A} >= 0, <= 1; # auxiliary variable for the linearisation

## objective function
minimize expectation: sum{p in 1..P, g in 1..G, a in A} y[p,g,a];

## constraints
subject to linearisation_1{p in 1..P, g in 1..G, a in A}:
	y[p,g,a] >= z[p,g,a] - sum{n in I_One[p,g,a]} x[n];
	
subject to linearisation_2{p in 1..P, g in 1..G, a in A, r in 1..T}:
	log(theta[r]) + 1/theta[r] * (z[p,g,a] - theta[r]) >= sum{n in I_OneHalf[p,g,a]} x[n] * log(0.5);
	
subject to population_size_conservation:
	sum{n in 1..N} x[n] = 2*N;
	
subject to population_size_equal_by_gender:
	sum{n1 in 1..Nm} x[n1] = sum{n2 in Nm+1..N} x[n2];
	
subject to x_upper_bound{n in 1..N}:
	x[n] <= x_UB;