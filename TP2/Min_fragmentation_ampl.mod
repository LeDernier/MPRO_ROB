param m;
param n;
param B;
param Amin;
param Amax;

param lambda;
param t;

set N := 1..n;
set M := 1..m;
set S := N cross M;
param C{S};

param d{(i1,j1) in S, (i2,j2) in S} := sqrt((i1-i2)^2 + (j1-j2)^2);

#Decision variables
var x{S} binary;
var y{S, S} binary;

minimize mean_distance:
	sum{(i1,j1) in S, (i2,j2) in S} (y[i1,j1,i2,j2] * d[i1,j1,i2,j2]) - (lambda * sum{(i,j) in S} (x[i,j]));

subject to surface_constraint_max:
	sum{(i,j) in S}(x[i,j]) <= Amax;

subject to surface_constraint_min:
	sum{(i,j) in S}(x[i,j]) >= Amin;
	
subject to cost_constraint:
	sum{(i,j) in S}(x[i,j] * C[i,j] * 10) <= B;

subject to neighbor_constraint{(i1,j1) in S, (i2,j2) in S}:
	x[i2,j2] >= y[i1,j1,i2,j2];
	
subject to one_distance_constraint{(i,j) in S}:
	sum{(i2,j2) in S: i2 != i or j2 != j}(y[i,j,i2,j2]) == x[i,j];
	


