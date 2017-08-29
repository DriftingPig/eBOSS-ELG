from math import *
def rom(a, b,f, eps = 1e-8):
    """Approximate the definite integral of f from a to b by Romberg's method.
    eps is the desired accuracy."""
    #   print 'a='+str(a)+' b='+str(b)
    #print ' f(a)='+str(f(a))+' f(b)='+str(f(b))
    R = [[0.5 * (b - a) * (f(a) + f(b))]]  # R[0][0]
    #print_row(R[0])
    n = 1
    while True:
        h = float(b-a)/2**n
        R.append((n+1)*[None])  # Add an empty row.
        
        R[n][0] = 0.5*R[n-1][0] + h*sum(f(a+(2*k-1)*h) for k in range(1, 2**(n-1)+1))  # for proper limits
        for m in range(1, n+1):
            R[n][m] = R[n][m-1] + (R[n][m-1] - R[n-1][m-1]) / (4**m - 1)
        #print_row(R[n])
        if abs(R[n][n-1] - R[n][n]) < eps:
            return R[n][n]
        n += 1

def dc(z): # Comoving distance from now to z
    c=2997.9
    h0=1
    return ((c/h0)*rom(0,z,evolution))

def evolution(z):
    om = float(0.31)
    ol=float(0.69)
    w=-1.
    fa=-3.*(1.+w)
    return 1./(sqrt(om*(1.+z)**3.+ol*(1.+z)**fa))

