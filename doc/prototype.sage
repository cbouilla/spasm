# -*- encoding: utf-8 -*-
K = GF(7)
n = 10 # lignes
m = 15 # colonnes
r = min(n, m) # borne sup sur le rang

# astuce pour générer une matrice avec un défaut de rang
A = random_matrix(K, m, n, density=0.33, sparse=True)
A[2] = A[1] + A[0]
A = A.T
A[3] = A[1] + A[2]

# initialise
p =  [ i for i in range(n) ] # permutation des lignes de L
qinv = [ -1 for i in range(m) ] # permutation des lignes de U
L = matrix(K, n, r)
U = matrix(K, r, m)
x = vector(K, [ 0 for i in range(m) ])
defficiency = 0

for i in range(n):
    # résoud x*U = A[i]
    a = copy( A[i] )
    for j in range(m):
        if a[j] != 0:
            if qinv[j] == -1:
                # pas de pivot colonne j --- identité implicite
                x[j] = a[j]
            else:
                x[j] = a[j] / U[ qinv[j], j ]
                a -= x[j] * U[ qinv[j] ]

    # recherche du pivot et enregistrement des coeffs de L
    ipiv = -1
    for j in range(m):
        if x[j] != 0:
            if qinv[j] == -1:
                if ipiv == -1:
                    ipiv = j
            else:
                L[i, qinv[j]] = x[j]

    if ipiv == -1:
        # pas de pivot sur cette ligne
        defficiency += 1
        p[n - defficiency] = i
    else:
        # pivot présent : l'enregistrer dans U et dans qinv
        assert( x[ipiv] != 0 )
        qinv[ipiv] = i - defficiency
        p[i - defficiency] = i
        U[i - defficiency, ipiv] = x[ipiv]
        L[i, i - defficiency] = 1

    # enregistrement des autres coeffs de U
    for j in range(m):
        if x[j] != 0:
            if qinv[j] == -1:
                U[i - defficiency, j] = x[j]
            x[j] = 0

    # on vérifie l'invariant
    assert (L*U)[:i] == A[:i]


# At this point, L*U == A
assert L*U == A

# now permute L and U to get a nice PLUQ decomposition of A

# P*L is lower-trapezoidal with unit "diagonal"
P = matrix(K, n, n, {(i, p[i]): 1 for i in range(n)}, sparse=True)

# U*Q.T is upper-trapezoidal with non-zero "diagonal"
Q = matrix(K, m, m, sparse=True)
k = 1
for i in range(m):
    if qinv[i] != -1:
        Q[qinv[i], i] = 1
    else:
        Q[m-k,i] = 1
        k += 1

LL = (P*L)[:,:n-defficiency] # keep only the first (n-defficiency) columns
UU = (U*Q.T)[:n-defficiency] # keep only the first (n-defficiency) rows

assert P.T*LL*UU*Q == A



def NNZ(A):
    return len( A.nonzero_positions() )

print "|A| = ", NNZ(A)
print "remplissage = ", (NNZ(L) + NNZ(U)) -  NNZ(A)
print "rank = ", n - defficiency
