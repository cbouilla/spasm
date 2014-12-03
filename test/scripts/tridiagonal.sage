#
# sage utility script to generate random sparse lower-triangular matrix
#
F = GF(2^16 + 1)
n = 100000

# lower-triangular
M = matrix(F, n, n, sparse=True)
M[0,0  ] = 1
M[0,1] = F.random_element()
for i in range(1, n-1):
    M[i,i-1] = F.random_element()
    M[i,i  ] = F.random_element()
    M[i,i+1] = F.random_element()
M[n-1,n-2] = F.random_element()
M[n-1,n-1] = 1

out = open("Matrix/tridiagonal", "w")
for (i,j) in M.nonzero_positions():
     out.write("{0} {1} {2}\n".format(i, j, M[i,j]))
out.close()
