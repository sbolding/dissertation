from scipy.special import expn
from scipy.integrate import quadrature


#integrate over a cell
x_start = 0.
x_end = 0.1/4.
sigma = 1000
psi_inc = 100000
q  = 1.0

result = quadrature(lambda x: psi_inc/(x_end-x_start)*expn(2,sigma*x),x_start,x_end,tol=1.e-12,maxiter=500)
print("HI")
psi_x1 = lambda x: (psi_inc*expn(2,sigma*x) + quadrature(lambda y:
    expn(1,sigma*y)*0.5*q,0,x,tol=1.e-12,maxiter=500)[0])
psi_x2 = lambda x: psi_inc*expn(2,sigma*x) + 0.5*q/sigma*(1 - expn(2,sigma*x))

print(psi_x1(0.001),psi_x1(0.01))
print(psi_x2(0.001),psi_x2(0.01))
