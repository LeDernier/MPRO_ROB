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

var x{NM} binary;
var y{(i,j) in NM,(k,l) in A[i,j]} binary;

# objective function

maximize obj: w1 * sum{(i,j) in NM} t[i,j]*(1-x[i,j]) + 
	w2*g*L*(sum{(i,j) in NM}(sum{(k,l) in A[i,j]} (x[i,j]-y[i,j,k,l]) 
	+ (4-card(A[i,j]))*x[i,j]));

subject to linearization{(i,j) in NM, (k,l) in A[i,j]}:
	y[i,j,k,l] >= x[i,j] + x[k,l] - 1;

	