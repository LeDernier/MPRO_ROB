param P;	 # ??
param m;
param n;
param w1;
param w2;
param L;
param g;

set N := 1..n;
set M := 1..m;
set NM := N cross M;
set A{i in N,j in M} := {(i2,j2) in NM: abs(i2-i) + abs(j2-j) = 1}; 			# neighbours of the zone (i,j)

param t{N,M};

# decision variables

var x{N,M} binary;
var d{N,M} >= 0;

# objective function

maximize obj: w1 * sum{(i,j) in NM} t[i,j]*(1-x[i,j]) + w2*g*L*sum{(i,j) in NM} (4*x[i,j] - d[i,j]);

subject to neighbours{(i,j) in NM}:
	d[i,j] >= sum{(k,l) in A[i,j]} x[k,l] - card(A[i,j])*(1-x[i,j]);
	
