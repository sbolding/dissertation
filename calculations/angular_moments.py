from scipy.special import expn
from scipy.integrate import quadrature


#integrate over a cell, do in terms of mfp
tau_start = 0.
tau_end = 1.0
sigma = 100
x_start = tau_start/sigma
x_end = tau_end/sigma
dx = x_end - x_start
psi_inc = 500
print("In MFP: ", sigma*(x_end - x_start))


#Compute values for multiple cells
for i in range(4):

    print(x_start,x_end)

    psi_x2 = lambda x: 6.*psi_inc*(expn(3,sigma*x) - 0.5*expn(2,sigma*x))
    phi_avg =psi_inc/(dx)*quadrature(lambda x: expn(2,sigma*x),x_start,x_end,tol=1.E-12,maxiter=500)[0]
    phi_eval = lambda x: 1/(sigma*(dx))*6.*psi_inc*(0.5*expn(3,sigma*x) - expn(4,sigma*x))
    phi_mu = phi_eval(x_end) - phi_eval(x_start)
    phi_x_int = lambda x: 6.*psi_inc/(dx*dx)*(x-0.5*(x_start+x_end))*expn(2,sigma*x)
    phi_x = quadrature(lambda x: phi_x_int(x),x_start,x_end,tol=1.E-12,maxiter=500)[0]
 
    print("---Results for %f < x < %f ---")
    print("Moments: ",phi_avg,phi_x,phi_mu)
    print("Corner values = ",phi_avg+phi_x-phi_mu,phi_avg-phi_x-phi_mu)
    print("Phi_avg outflow = ",phi_avg+phi_x)
    print("")

    x_start += dx
    x_end += dx



