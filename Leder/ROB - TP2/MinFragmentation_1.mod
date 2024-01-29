param m;
param n;		# hypothesis: n = m
param B;
param Amin;
param Amax;

set N := 1..n;
set M := 1..m;
set P := N cross M;

param c{N,M};
param lambda;

param distance{(i1,j1) in P, (i2,j2) in P} := abs(i2 - i1) + abs(j2 - j1);

# decision variables

var x{P} binary;
var y{P,P} binary;

# objective function

minimize obj: sum{(i1,j1) in P, (i2,j2) in P: i2 != i1 or j2 != j1} distance[i1,j1,i2,j2]*y[i1,j1,i2,j2] - lambda*sum{(i,j) in P} x[i,j];

# constraints
subject to PPV_1{(i1,j1) in P}:
	sum{(i2,j2) in P: i2 != i1 or j2 != j1} y[i1,j1,i2,j2] = x[i1,j1];
	
subject to PPV_2{(i1,j1) in P, (i2,j2) in P: i2 != i1 or j2 != j1}:
	x[i2,j2] >= y[i1,j1,i2,j2];
	
subject to area_min:
	sum{(i,j) in P} x[i,j] >= Amin;			# hypothesis: a[(i,j)] = 1 km
	
subject to area_max:
	sum{(i,j) in P} x[i,j] <= Amax;			# hypothesis: a[(i,j)] = 1 km
	
subject to cost_max:
	sum{(i,j) in P} c[i,j]*x[i,j] <= B;

