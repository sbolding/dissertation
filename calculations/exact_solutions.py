from scipy.special import expn
from scipy.integrate import quadrature


#integrate over a cell
x_start = 0.
x_end = 1./20.
sigma = 1000
psi_inc = 50
q  = 1.0

psi_x1 = lambda x: (psi_inc*expn(2,sigma*x) + quadrature(lambda y:
    expn(1,sigma*y)*0.5*q,0,x,tol=1.e-12,maxiter=500)[0])
psi_x2 = lambda x: psi_inc*expn(2,sigma*x) + 0.5*q/sigma*(1 - expn(2,sigma*x)) + 0.5*q/sigma
phi_avg = 1/(x_end-x_start)*quadrature(lambda x: psi_x2(x),x_start,x_end,tol=1.E-12,maxiter=500)[0]

print(phi_avg, psi_x2(0.5*(x_end+x_start)),psi_x2(x_end))
