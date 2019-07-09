from __future__ import print_function
from dolfin import *
from fenics import *
import matplotlib.pyplot as plt
from mshr import *
import numpy as np
import random
import csv


L,H,r = 400.0,400.0,20.0
k,w,M,theta = 0.29,0.1,5,0.5
a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10 = 0,0,8.072789087,-81.24549382,408.0297321,-1244.129167,2444.046270,-3120.635139,2506.663551,-1151.003178,230.2006355
#T = 50
#num_steps = 100
dt = 0.1
mesh = RectangleMesh(Point(-L/2,-H/2),Point(L/2,H/2),200,200,"crossed")
facet = MeshFunction("bool", mesh, 0)
facet.set_all(False)

class Inside(SubDomain):
       def inside(self,x,on_boundary):
              return near(sqrt(pow((x[0]-0),2)+pow((x[1]-0),2)),20,10)#sqrt(pow((x[0]-400/2),2)+pow((x[1]-400/2),2)) <=(20+10)
Inside().mark(facet, True)
mesh = refine(mesh, facet)


class CahnHilliard(NonlinearProblem):
    def __init__(self,a,L):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
    def F(self,b,x):
        assemble(self.L,tensor=b)
        #[bc.apply(b) for bc in bcs]

    def J(self,A,x):
        assemble(self.a,tensor=A)
       # [bc.apply(A) for bc in bcs]
tol1 = 1e-14
a = r/0.9
b = r*0.9
eta0 = Expression(("sqrt(pow((x[0]-0),2)/pow(a,2)+pow((x[1]-0),2)/pow(b,2)) <= (1+tol1) ? 1 : 0.0065","0"),degree =0,L=L,H=H,r=r,tol1 = tol1,a=a,b=b)

#eta0 = Expression(("sqrt(pow((x[0]-L/2),2)+pow((x[1]-H/2),2)) <= r ? 1 :0.0065","0"),degree =0,L=L,H=H,r=r)
d =1
Vne = FiniteElement('CG',mesh.ufl_cell(),d) #for eta 
Vme = FiniteElement('CG',mesh.ufl_cell(),d) #for mu
V = FunctionSpace(mesh,MixedElement([Vne,Vme]))
V1,V2 = TestFunction(V)
VTrial= TrialFunction(V)
(deta,dmu) = split(VTrial)
VOld = Function(V)

VNew = Function(V)
eta_init = eta0
VOld.interpolate(eta_init)
(eta_0,mu_0) = split(VOld)
VNew.interpolate(eta_init)
(eta,mu) = split(VNew)
##compute dfdn ###
eta = variable(eta)
f1 = w*(a0 +a1*eta+a2*eta**2+a3*eta**3+a4*eta**4+a5*eta**5+a6*eta**6+a7*eta**7+a8*eta**8+a9*eta**9+a10*eta**10)
dfdn = diff(f1, eta)

##compute elasto energy ####
def elasticE (eta):  
      n = variable(eta)
      C11 = 250.
      C12 = 150.
      C66 = 100.
      Cm = as_matrix([[C11,C12,0],[C12,C11,0],[0,0,C66]])
      Cp = 1*Cm
      h = n**3*(6*n**2 -15*n + 10)
      C = Cm*(1-h) + Cp*h
      eigen  = 0.005*h
      class xaxis(SubDomain):
          def inside(self,x,on_boundary):
              return (near(x[0],0.) and near (x[1],-H/2)) or (near(x[0],0.) and near (x[1],0.)) or (near(x[0],0.) and near (x[1],H/2))

      class yaxis(SubDomain):
          def inside(self,x,on_boundary):
              return (near(x[0],-L/2) and near (x[1],0.)) or (near(x[0],0.) and near (x[1],0.)) or (near(x[0],L/2) and near (x[1],0.))
      facets = MeshFunction("size_t", mesh, 1)
      facets.set_all(0)
      xaxis().mark(facets, 1)
      yaxis().mark(facets, 2)
      Vu = VectorFunctionSpace(mesh,'Lagrange',2) 
      bc1 = DirichletBC(Vu.sub(1), Constant(0),facets,1)
      bc2 = DirichletBC(Vu.sub(0), Constant(0),facets,2)
      bcs = [bc1,bc2]
      f = Constant((0.,0.))
      def eps(v):
          return sym(grad(v))
      def strain2voigt(e):
          return as_vector([e[0,0],e[1,1],2*e[0,1]])
      def voigt2stress(s):
          return as_tensor([[s[0],s[2]],[s[2],s[1]]])
      def est(eigen):
          return as_vector([eigen,eigen,0])    
      def sigma0(eigen):
          return voigt2stress(dot(C,est(eigen)))
      def sigma(v):
          return voigt2stress(dot(C,strain2voigt(eps(v))))
      def sigmar(v,eigen):
          return (sigma(v) - sigma0(eigen))   
      def stress2voigt(p):
          return as_vector([p[0,0],p[1,1],p[0,1]])
      def el(v,eigen):
          return inner(sigmar(v,eigen),sym(grad(v)))*dx
      du = TrialFunction(Vu)
      dv = TestFunction(Vu)
      Wint = inner(sigmar(du,eigen),grad(dv))*dx
      aM = lhs(Wint)
      LM = rhs(Wint) + inner(f, dv)*dx
      u = Function(Vu,name="Displacement")
      solve(aM == LM,u,bcs,solver_parameters={'linear_solver':'mumps'})
      u = project(u,Vu)
      strain = sym(grad(u)) - voigt2stress(est(eigen))
      E = 0.5*inner(sigmar(u,eigen),strain)
      print(assemble(E*dx))
      return E     
solver = NewtonSolver()
solver.parameters["linear_solver"]="lu"
solver.parameters["convergence_criterion"]="incremental"
solver.parameters["relative_tolerance"]= 1e-6
file1 = File ("sol_4e4_p01_em7/solution_.pvd","compressed")
t = 0.0
file1 << (VNew.split()[0], t)
#mu_0 = dfdn + elasticE(eta) - k*(div(grad(eta)))
while True:
    t += dt
    VOld.vector()[:] = VNew.vector()
    El = assemble(elasticE(eta)*dx)
    Bl = assemble(f1*dx)
    Total_Energy = (El+Bl+assemble((k/2*dot(grad(eta),grad(eta)))*dx))
    print(Total_Energy)
    TE = [t,El,Bl,Total_Energy-El-Bl,Total_Energy]

    mu_mid = (1.0-theta)*mu_0 + theta*mu
    cahn1 = (eta - eta_0)/dt*V1*dx + M*dot(grad(mu_mid),grad(V1))*dx
    cahn2 = mu*V2*dx -diff(elasticE(eta),(eta))*V2*dx-dfdn*V2*dx -k*dot(grad(eta),grad(V2))*dx
    cahn = cahn1+cahn2
    a = derivative(cahn,VNew,VTrial)
    problem = CahnHilliard(a,cahn)
    (no_of_iterations,converged) = solver.solve(problem,VNew.vector())
    if no_of_iterations <4:
        dt = 2*dt
    else:
        dt = dt/2
    file1 << (VNew.split()[0], t)
    q = assemble(abs(mu - mu_0)*dx)/len(mesh.coordinates())
    with open('EnergyN_dt_p001.csv',mode='a') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(TE)
    print(q)
    if (q < 1e-10):
        break

