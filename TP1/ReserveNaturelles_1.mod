param n;
param m;
param p;
param q;

set N := 1..n;			# horizontal index (grid)
set M := 1..m;			# vertical index (grid)
set P := 1..p;			# in-danger speces
set Q := p+1..p+q;		# common speces
set PQ := P union Q;
set V{i in N,j in M} := {(i2,j2) in N cross M: (abs(i2-i) = 1 and abs(j2 - j) <= 1) or (abs(i2-i) <= 1 and abs(j2 - j) = 1)}; 			# neighbours of the zone (i,j)

param c{N,M};					# costs
param proba{PQ,N,M} default 0;	# survival probability in a zone if the zone is protected
param alpha{PQ};				# survival probability threshold 

# decision variables

var x{N,M} binary default 0;
var y{N,M} binary default 0;



# objective function

minimize f: sum{i in N, j in M}  c[i,j]*x[i,j];

subject to central_zone{i1 in N, j1 in M}:
	sum{(i2,j2) in V[i1,j1] union {(i1,j1)}} x[i2,j2] >= card(V[i1,j1] union {(i1,j1)})*y[i1,j1];	# implies y[i,j] = 1 if the surrounding zones are selected (x|i',j'] = 1), even if x[i,j] = 0

/* Remark: this constraint is contained in the one above (we replaced "V[i1,j1]" by "V[i1,j1] union {(i1,j1)}}")
subject to chosen_zone{i in N, j in M}:	# implies x[i,j] = 1 whenever y[i,j] = 1
	x[i,j] >= y[i,j]; 
*/

subject to prob_extinction_esp_in_danger{k in P}:
	sum{i in N, j in M} log(1-proba[k,i,j])*y[i,j] <= log(1-alpha[k]);

subject to prob_extinction_esp_commons{k in Q}:
	sum{i in N, j in M} log(1-proba[k,i,j])*x[i,j] <= log(1-alpha[k]);
