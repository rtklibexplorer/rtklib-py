"""
integer ambiguity resolution by LAMBDA

reference :
     [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
         a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82,
         1995
     [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
         integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
         
Copyright (c) 2021 Rui Hirokawa (from CSSRLIB)
Copyright (c) 2022 Tim Everett

"""

import numpy as np
from numpy.linalg import inv


def LD(Q):
    """ LD factorization (Q=L'*diag(D)*L) """
    n = len(Q)
    L = np.zeros((n, n))
    d = np.zeros(n)
    A = Q.copy()
    for i in range(n-1, -1, -1):
        d[i] = A[i,i]
        if d[i] <= 0.0:
            print('LD Factorization error')
            raise SystemExit
        L[i,:i+1] = A[i,:i+1] / np.sqrt(d[i])
        for j in range(i):
            A[j,:j+1] -= L[i,:j+1] * L[i,j]
        L[i,:i+1] /= L[i,i]

    return L, d


def reduction(L, d):
    """ lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) """
    n = len(d)
    Z = np.eye(n)
    j = k = n-2
    while j >= 0:
        if j <= k:
            for i in range(j+1, n):
                # L,Z=gauss(L,Z,i,j)
                mu = round(L[i,j])
                L[i:,j] -= mu * L[i:,i]
                Z[:,j] -= mu * Z[:,i]

        delta = d[j] + L[j+1,j]**2 * d[j+1]
        if delta + 1e-6 < d[j+1]:  # compared considering numerical error 
            eta = d[j] / delta
            lam = d[j+1] * L[j+1,j] / delta
            d[j] = eta * d[j+1]
            d[j+1] = delta
            L[j:j+2,:j] = np.array([[-L[j+1, j],1], [eta,lam]]) @ L[j:j+2,:j]
            L[j+1,j] = lam
            # swap L's and Z's
            L[j+2:,j], L[j+2:,j+1] =  L[j+2:,j+1].copy() , L[j+2:,j].copy()
            Z[:,j], Z[:,j+1] = Z[:,j+1].copy(), Z[:,j].copy()

            j, k = n -2, j

        else:
            j -= 1
    return L, d, Z


def search(L, d, zs, m=2):
    """ modified lambda (mlambda) search (ref. [2])
* args     m      I  number of fixed solution
           L,d    I  transformed covariance matrix
           zs     I  transformed double-diff phase biases
           zn     O  fixed solutions
           s      O  sum of residuals for fixed solutions """

    n = len(d)
    nn = 0
    imax = 0
    Chi2 = 1e18
    S = np.zeros((n, n))
    dist = np.zeros(n)
    zb = np.zeros(n)
    z = np.zeros(n)
    step = np.zeros(n)
    zn = np.zeros((n, m))
    s = np.zeros(m)
    k = n - 1
    zb[-1] = zs[-1]
    z[-1] = round(zb[-1])
    y = zb[-1] - z[-1]
    step[-1] = np.sign(y) # step towards closest integer
    if step[-1] == 0:
        step[-1] = 1
    for _ in range(10000):
        # newdist=sum(((z(j)-zb(j))^2/d(j)))
        newdist = dist[k] + y**2 / d[k]
        if newdist < Chi2:
            # Case 1: move down
            if k != 0:
                k -= 1
                dist[k] = newdist
                S[k, :k+1] = S[k+1, :k+1] + (z[k+1] - zb[k+1]) * L[k+1,:k+1]
                zb[k] = zs[k] + S[k,k]
                z[k] = round(zb[k])
                y = zb[k] - z[k]
                step[k] = np.sign(y) 
                if step[k] == 0:
                    step[k] = 1
            # Case 2: store the found candidate and try next valid integer
            else:
                if nn < m: # store the first m initial points
                    if nn == 0 or newdist > s[imax]:
                        imax = nn
                    zn[:,nn] = z
                    s[nn] = newdist
                    nn += 1
                else:
                    if newdist < s[imax]:
                        zn[:,imax] = z
                        s[imax] = newdist
                        imax = np.argmax(s)
                    Chi2 = s[imax]
                z[0] += step[0] # next valid integer
                y = zb[0]-z[0]
                step[0] = -step[0] - np.sign(step[0])
        # Case 3: exit or move up
        else:
            if k == n - 1:
                break            
            k += 1 # move up
            z[k] += step[k] # next valid integer
            y = zb[k] - z[k]
            step[k] = -step[k] - np.sign(step[k])

    order = np.argsort(s) # sort by s
    s = s[order]
    zn = zn[:, order]

    return zn, s


def mlambda(a, Q, m=2):
    """lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2])
* args     m      I  number of fixed solutions
*          a      I  float parameters (n x 1) (double-diff phase biases)
*          Q      I  covariance matrix of float parameters (n x n)
*          afix_  O  fixed solutions (n x m)
*          s      O  sum of squared residulas of fixed solutions (1 x m) """

    # LD (lower diaganol) factorization (Q=L'*diag(D)*L) 
    L, d = LD(Q)
    L, d, Z = reduction(L, d)
    invZt = np.round(inv(Z.T))
    z = Z.T @ a
    # mlambda search 
    #        z = transformed double-diff phase biases
    #        L,D = transformed covariance matrix
    E, s = search(L, d, z, m)
    afix_ = invZt @ E
    return afix_, s
