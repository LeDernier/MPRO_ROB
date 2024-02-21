param P;
param m;
param n;
param w1;
param w2;
param L;
param g;

set M := 1..m;
set N := 1..n;
set S := N cross M;
set A{(i,j) in S} := {(i2,j2) in S: abs(i2-i) + abs(j2-j) = 1};

param t{S};

#decision variables
var x{S} >= 0, <= 1;
var y{S,S} >= 0, <= 1;

maximize p1:
	w1*sum{(i,j) in S}(t[i,j]*(1-x[i,j])) 
	+ w2*g*L*sum{(i,j) in S}
	(sum{(k,l) in A[i,j]}(x[i,j] - y[i,j,k,l])
	+(4-card(A[i,j]))*x[i,j]);
	
subject to calc_y{(i,j) in S,(k,l) in A[i,j]}:
	y[i,j,k,l] >= x[i,j] + x[k,l] - 1;

subject to min_60:
	sum{(i,j) in S}(x[i,j]) >= 60;