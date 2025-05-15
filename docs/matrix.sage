# This script provides checks of the computations in.
# "A projective variety with discrete, non-finitely generated automorphism group", J. Lesieutre

# The file has been checked to run in SageMath 7.5.1.
# It can be run from the command line via "sage final.sage" and will print the most important information to stdout.
# Many auxiliary computations are stored in variables and can be accessed by running the script in interactive mode.

# The output of the script is attached to the paper in the file "final.sage.out".

########################

# First name the vectors in a basis for N^1(S) to match the notations of the paper.
(H, E01, E02, E03, E04, E05, E12, E13, E14, E15, E23, E24, E25, E34, E35, E45) = matrix.identity(16)

# Now define the classes of the 6 lines
(L0, L1, L2, L3, L4, L5) = (H - E01 - E02 - E03 - E04 - E05,
                            H - E01 - E12 - E13 - E14 - E15,
			    H - E02 - E12 - E23 - E24 - E25, 
                            H - E03 - E13 - E23 - E34 - E35,
			    H - E04 - E14 - E24 - E34 - E45, 
                            H - E05 - E15 - E25 - E35 - E45)
L = matrix([L0,L1,L2,L3,L4,L5])			    
# and the anticanonical class.
K = (3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1)

# Define the intersection form on N^1(S)
intersect_N1 = diagonal_matrix([1]+[-1]*15)

# It will be useful to get the class Eij from i and j, which is done by the following functions.
# For example, exceptional_index([1,2]) = 6 since E12 is the 6th entry.
# also define Ex([1,2]), which returns the class E12.
def exceptional_index(v):
  i = v[0]
  j = v[1]
  if(i > j):
    temp = i
    i = j
    j = temp
  return j-(i*(i-9))/2

def Ex(v):
  i = v[0]
  j = v[1]
  if(i != j):
    return matrix.identity(16).column(exceptional_index(v))

# We also define function to get the class of the line Li.
# L1 is position 0 in the list L, so decrease index by 1.
def Li(i):
  return L[i]

# The first step is to compute the matrix for the action of the map s_0 on N^1(S)

# There are many classes whose image under s_0 is easy to see geometrically, and
# these provide enough information to determine the matrix.

# The following function takes as input the indices of the six lines being used.
# [La,Lb,Lc] is the first triple, with La distinguished.
# [Ld,Le,Lf] is the first triple, with Ld distinguished.
# So the map used in the paper is [a,b,c,d,e,f] = [0,1,4,5,2,3] (or [0,2,3,5,1,4])
 
def get_matrix_aux(v):

  [a,b,c,d,e,f] = v

  # The list s0_Images contains entries of the form [D,s0(D)] for which the
  # class s0(D) can be seen geometrically.
  # Since it's an involution, it goes the other way too.
  s0_Images = ( [ Li(a), Li(a) ], [ Li(d), Li(d) ],
                [ Li(b), Li(c) ], [ Li(e), Li(f) ],
                [ Ex([a,b]), Ex([a,c]) ], [ Ex([a,e]), Ex([a,f]) ],
		[ Ex([d,b]), Ex([d,c]) ], [ Ex([d,e]), Ex([d,f]) ],
		[ Ex([b,c]), Ex([b,c]) ], [ Ex([e,f]), Ex([e,f]) ],
                [ Ex([b,e]), H-Ex([a,d])-Ex([b,e]) ],
		[ Ex([b,f]), H-Ex([a,d])-Ex([b,f]) ],
		[ Ex([c,e]), H-Ex([a,d])-Ex([c,e]) ],
		[ Ex([c,f]), H-Ex([a,d])-Ex([c,f]) ],
		[ K, K ] )
  # Now pull off the first and second entries of the list.		
  entry_1 = matrix([ s0_Images[j][0] for j in range(len(s0_Images))] ).transpose()
  entry_2 = matrix([ s0_Images[j][1] for j in range(len(s0_Images))] ).transpose()
  
  # Let X be a matrix whose columns are the classes are being mapped, and B what they are mapped to.
  # entry_1 maps to entry_2 and vice-versa, so we want MX = B where:
  X = block_matrix(1,2,[entry_1,entry_2])
  B = block_matrix(1,2,[entry_2,entry_1])
  
  # To arrange that MX = B, take M = (BX^T)*(XX^T)^{-1} (which is defined, since X has full rank)
  M = (B*X.transpose())*(X*X.transpose()).inverse()

  return (M,X,B)

def get_matrix(v):
# Same function, but don't return known classes X and B, only the matrix M.

  return get_matrix_aux(v)[0]



# Now set which lines we are actually using and generate the corresponding matrix.
m1 = [0,1,4,5,2,3]
(M,X,B) = get_matrix_aux(m1)

print("Labelling the lines as "+"(*L"+str(m1[0])+",L"+str(m1[1])+",L"+str(m1[2])+"),"+
                                "(*L"+str(m1[3])+",L"+str(m1[4])+",L"+str(m1[5])+")"+".\n")

print("Matrix for action of s_0*:")
print(M)
print("")

# We now carry out some checks to make sure that the matrix M behaves as expected.
print("Check for correct action on known classes:  M*X == B:")
print(M*X == B)
print("Check that M*M == id:")
print(M*M == matrix.identity(16))
print("Check that M preserves intersection form:")
print(M.transpose()*intersect_N1*M == intersect_N1)
print("")

# We are interested in the classes which are:
# (1) Contained in the 1-eigenspace
# (2) Have the same intersections with L_i as the class H-E12.
# (3) Positive intersection with known irreducible effective classes (or in the list)
# (4) D^2 = 0

# We first construct a matrix containing a basis for the 1-eigenspace of M as its rows.
eigenspace_plus = (M-1*matrix.identity(16)).right_kernel().basis_matrix()

# We will now compute the polyhedron containing all classes which are contained in
# the 1-eigenspace and have suitable intersections with the lines L_i.
# To work with polyhedra it will be useful to characterize the 1-eigenspace by a set of equalities
# on the coefficients.  Our class D is supposed to be contained in the 1-eigenspace,
# which means it's orthogonal to each element of the kernel of eigenspace_plus.

# To that end, we generate a matrix whose rows are orthogonal to the 1-eigenspace.
# We then add a column of 0's to the left of this matrix for use with Sage's Polyhedron object.
eigenspace_plus_complement = eigenspace_plus.right_kernel().basis_matrix()
eigenspacezeros = [0]*eigenspace_plus_complement.nrows()
subspace_eqnsT = eigenspace_plus_complement.columns()
subspace_eqnsT.insert(0,eigenspacezeros)
subspace_eqns = matrix(subspace_eqnsT).transpose()

# Now, the equalities that force a divisor to have the desired intersections.
# If the marked lines are L0 and L5, the preserved fibration is H-E05.
# This has intersection 0 with L0 and L5, and 1 with the others.
# Should be [0,1,1,1,1,0] if L0 and L5 are marked.
nefclass = H - Ex([m1[0],m1[3]])
Li_intersections = L*intersect_N1*nefclass

# Now set up the equations to have D.Li as desired.
Li_eqnsT = (L*intersect_N1).columns()
Li_eqnsT.insert(0,-1*vector(Li_intersections)) # need the - sign because Sage writes Ax+b=0 instead of Ax=b.
Li_eqns = matrix(Li_eqnsT).transpose()

# Combine the equalities into a single list.
intersection_eqns = block_matrix([subspace_eqns,Li_eqns],nrows=2,ncols=1)

# Now write a list of known effective classes, as a list.
# These give us the inequalities: our classes must have positive intersection with these.
# The first of these are the 15 (-1)-curves created by blowing up
exceptionals = matrix.identity(16)[:,1:16].columns()
# Now make a list of the two (-1)-curves we care about
# For the labelling m1, these are the two classes H-E05-E23 and H-E05-E14.
lines = [ H - Ex([m1[0],m1[3]]) - Ex([m1[1],m1[2]]) , H - Ex([m1[0],m1[3]]) - Ex([m1[4],m1[5]]) ]

# Finally, form a matrix effectives whose rows are all the known effective classes.
effectives = exceptionals + lines

# Now, convert this into a matrix encoding the inequalities for a divisor D to have
# positive intersection with all these effective classes.
effectives_list = matrix(effectives)
effectives_ineqsT = (effectives_list*intersect_N1).columns()
effectives_ineqsT.insert(0,[0]*effectives_list.nrows())
effectives_ineqs = matrix(effectives_ineqsT).transpose()

# Generate the polyhedron of points satisfying all the constraints.
myregion = Polyhedron(ieqs=effectives_ineqs,eqns=intersection_eqns)
myvertices = matrix(myregion.vertices_list()).transpose()
myrays = matrix(myregion.rays_list()).transpose()

generating_classes = block_matrix([[myvertices,myrays]])
generator_intersections = generating_classes.transpose()*intersect_N1*generating_classes

print("Vertices of polyhedron:")
print(str(myvertices.transpose()))
print("Rays of polyhedron:")
print(str(myrays.transpose()))
print("")
print("Intersection numbers of generating classes are:")
print(generator_intersections)
print("\n")


# Now look at a map that fixes the line pointwise.
# We'll use the commutator of mu(0) and mu(1), but to do this, we first need matrices for sigma and tau.
M1325 = get_matrix([0,1,3,4,2,5])
M1345 = get_matrix([0,1,3,2,4,5])
M1234 = get_matrix([0,1,2,5,3,4])
Ms = M1234*M1345*M1325 # Matrix for sigma

M1523 = get_matrix([0,1,5,4,2,3])
M2435 = get_matrix([0,2,4,1,3,5])
Mt = M1345*M2435*M1523 # Matrix for tau

def mu(n):
  return (Mt^(-n))*Ms*(Mt^n)

def commutator(A,B):
  return A*B*(A.inverse())*(B.inverse())

def specrad(A):
  return max( [abs(x) for x in A.eigenvalues()] )

compose = commutator(mu(0),mu(1))
degree = compose[0][0]
lambda1 = max( [abs(x) for x in compose.eigenvalues()] )

print("Matrix for transformation fixing L_0 pointwise")
print("Degree of transformation:")
print(degree)
print("Dynamical degree of transformation:")
print(lambda1)