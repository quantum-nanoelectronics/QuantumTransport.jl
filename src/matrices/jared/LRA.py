from scipy.sparse import diags
import numpy as np
from scipy.sparse.linalg import eigs

offset = lambda arr: np.complex64(arr) + complex(0,1e-5)

def truncate_eigens(eigen_vals, eigen_vecs, cutoff):
    for i,e in enumerate(eigen_vals):
        if ( cmplx_mag(e) > cutoff ):
            eigen_vecs[i] = np.zeros( len(np.transpose(eigen_vecs)) )
            # remove eigen energy too?
            eigen_vals[i] = 0
    return eigen_vals, eigen_vecs

def construct_laplace(n):
    return (np.diag(offset(np.ones(n-1)) , -1) +
            np.diag(offset(np.full(n,-2)))     +
            np.diag(offset(np.ones(n-1)), 1) )

def construct_pos(n):
    return (np.diag( [np.sqrt(complex(i+1,0)) for i in range(n-1)], -1) +
            np.diag( [np.sqrt(complex(i+1,0)) for i in range(n-1)], 1) )

def cmplx_mag(u):
    return np.sqrt(u * np.conjugate(u))

n = 4
H = construct_laplace(n)
E, phi = truncate_eigens(*eigs(H), 3)
S = phi # matrix of eigen vectors as columns
X = construct_pos(n)
S_dagger = np.conjugate( np.transpose(S) )
Xs = S_dagger * X * S
xs , psi = eigs(Xs)

I = np.identity(n)
E_ieta = offset(E)
G_R = (E_ieta) * I - H # minus self energy terms
