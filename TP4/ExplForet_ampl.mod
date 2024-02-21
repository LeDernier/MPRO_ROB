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
var x{S} binary;
var d{S} >= 0;

maximize p1:
	w1*sum{(i,j) in S}(t[i,j]*(1-x[i,j])) + w2*g*L*sum{(i,j) in S}(4*x[i,j]-d[i,j]);
	
subject to calc_d{(i,j) in S}:
	d[i,j] >= sum{(k,l) in A[i,j]}(x[k,l]) - card(A[i,j])*(1-x[i,j]);

subject to min_60:
	sum{(i,j) in S}(x[i,j]) >= 60;
