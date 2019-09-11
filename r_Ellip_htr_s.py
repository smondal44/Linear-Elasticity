"""
Phase field community hub developed by Center for Hierarchial Materials Design (CHiMaD) and the National Institute of 
Standards and Technology (NIST) has published several problems along with standard solutions in order new implementations
of phase field models can be benchmarked. In this work, we present the FEM based FEniCS implementation to solve the 
Linear elasticity problem for elliptical heterogeneous precipitate of radious 20nm in the domain of 400*400. 
This work has been completed during summer'2019 at MATERIALS AND PROCESS MODELLING LAB
of IIT Bombay.
Author: Shubhajit Mondala(a), Sushil Kumar(b), M.P. Gururajan(b)
(a) Department of Minerals, Metallurgical and Materials Engineering, Indian Institute of Technology Bhubaneswar, Khurda 752050,India
(b) Department of Metallurgical Engineering and Materials Science,Indian Institute of Technology Bombay, Powai, Mumbai 400076,India
"""
from __future__ import print_function
from dolfin import *
from fenics import *
import matplotlib.pyplot as plt
from mshr import *
import numpy as np
import random
import csv

L,H,r = 400.0,400.0,20.0 # domain size and radious of the ppt.
k,w,M,theta = 0.29,0.1,5,0.5 #value of the constants used for model formulation
# coefficient of the polynomial of the bulk free energy
a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10 = 0,0,8.072789087,-81.24549382,408.0297321,-1244.129167,2444.046270,-3120.635139,2506.663551,-1151.003178,230.2006355 
dt = 0.1 #initial delta dt
mesh = RectangleMesh(Point(-L/2,-H/2),Point(L/2,H/2),200,200,"crossed") # generation of the rectangular mesh
# Refinig the mesh close to the peimiter of the circular ppt.
facet = MeshFunction("bool", mesh, 0)
facet.set_all(False)
class Inside(SubDomain):
       def inside(self,x,on_boundary):
              return near(sqrt(pow((x[0]-0),2)+pow((x[1]-0),2)),20,10)#sqrt(pow((x[0]-400/2),2)+pow((x[1]-400/2),2)) <=(20+10)
Inside().mark(facet, True)
mesh = refine(mesh, facet)
# the class is for assembling the jacobian and non-linear form of the problem
class CahnHilliard(NonlinearProblem):
    def __init__(self,a,L):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
    def F(self,b,x):
        assemble(self.L,tensor=b)
    def J(self,A,x):
        assemble(self.a,tensor=A)
tol1 = 1e-14 #tolerance to approximate circular region more accurately
# Expression for setting the initial condition of the conserved -order parameter.
a = r/0.9 #major axis of the ellipse
b = r*0.9 #minor axis of the ellipse
eta0 = Expression(("sqrt(pow((x[0]-0),2)/pow(a,2)+pow((x[1]-0),2)/pow(b,2)) <= (1+tol1) ? 1 : 0.0065","0"),degree =0,L=L,H=H,r=r,tol1 = tol1,a=a,b=b)
d =1 #degree of the function-space
Vne = FiniteElement('CG',mesh.ufl_cell(),d) #finite element for scalar eta 
Vme = FiniteElement('CG',mesh.ufl_cell(),d) #finite element for scalar mu
V = FunctionSpace(mesh,MixedElement([Vne,Vme])) # mixing the two scalar finite element to make mixed function spcae,here sequence is importanrt.
V1,V2 = TestFunction(V)
VTrial= TrialFunction(V)
(deta,dmu) = split(VTrial)
VOld = Function(V)
VNew = Function(V)
VOld.interpolate(eta0)
(eta_0,mu_0) = split(VOld)
VNew.interpolate(eta0)
(eta,mu) = split(VNew)
#Compute the variational derivative of the bulk-free energy
eta = variable(eta)
f1 = w*(a0 +a1*eta+a2*eta**2+a3*eta**3+a4*eta**4+a5*eta**5+a6*eta**6+a7*eta**7+a8*eta**8+a9*eta**9+a10*eta**10)
dfdn = diff(f1, eta)
#Compute elastic energy of the system by solving the equation of the mechanical equilibrium
def elasticE (eta):  
      n = variable(eta)
#the value of 3 independent elastic tensor as considering cubic anisotropy
      C11 = 250.0 
      C12 = 150.0
      C66 = 100.0
      Cm = as_matrix([[C11,C12,0],[C12,C11,0],[0,0,C66]]) # the elastic tensor for matrix in voigt notation
      Cp = 1.1*Cm #considering the heterogeneous ppt.
      h = n**3*(6*n**2 -15*n + 10) #wang polynomial to interpolate the variable smoothly from one phase to another phase across the smooth boundary
      C = Cm*(1-h) + Cp*h
      eigen  = 0.005*h
#the following two classes for applying the boundary condition in order to remove the null space of the solution
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
      Vu = VectorFunctionSpace(mesh,'Lagrange',2)  # formualtinf vectorfunction-space for solving the mechanical equilibrium equation
      bc1 = DirichletBC(Vu.sub(1), Constant(0),facets,1)
      bc2 = DirichletBC(Vu.sub(0), Constant(0),facets,2)
      bcs = [bc1,bc2]
      f = Constant((0.,0.)) # to make body force of the system is zero
      def eps(v): # return the strain in tensor form taking displacement as input
          return sym(grad(v))
      def strain2voigt(e): #converting the strain in tensor form to the strain in voigt notation
          return as_vector([e[0,0],e[1,1],2*e[0,1]])
      def voigt2stress(s): #Converting the stress in voigt notation to the tensor form
          return as_tensor([[s[0],s[2]],[s[2],s[1]]])
      def est(eigen): #incorporating the eigen strain in oigt notation format
          return as_vector([eigen,eigen,0])    
      def sigma0(eigen): #calculating the stress due to eigen strain
          return voigt2stress(dot(C,est(eigen)))
      def sigma(v): #calculating the stress due to constrained displacement
          return voigt2stress(dot(C,strain2voigt(eps(v))))
      def sigmar(v,eigen): #resultant stress = stress due to constrained displacement - stress due to eigen strain
          return (sigma(v) - sigma0(eigen))   
#defining trial and test function for solving the equation of mechanical equilibrium
      du = TrialFunction(Vu)
      dv = TestFunction(Vu)
      Wint = inner(sigmar(du,eigen),grad(dv))*dx #variational form of the mechanical equilibrium equation
      aM = lhs(Wint) # to separate the bi-linear form the variational form(Wint)
      LM = rhs(Wint) + inner(f, dv)*dx 
      u = Function(Vu,name="Displacement")
      solve(aM == LM,u,bcs,solver_parameters={'linear_solver':'mumps'})
      u = project(u,Vu) # to interpolate/project the solution to the function of function -space(Vu)
      strain = sym(grad(u)) - voigt2stress(est(eigen)) #elastic energy = total strain- eigen strain
      E = 0.5*inner(sigmar(u,eigen),strain)
      print(assemble(E*dx))  # priniting the elastic energy of the system at each time step
      return E     # returning the elastic energy
#Specifying the solver specification to have control over the solver
solver = NewtonSolver()
solver.parameters["linear_solver"]="lu"
solver.parameters["convergence_criterion"]="incremental"
solver.parameters["relative_tolerance"]= 1e-6
file1 = File ("sol_folder/solution_.pvd","compressed") # initializing the file name to store the solution
t = 0.0 #initial time
file1 << (VNew.split()[0], t) # initial value of the conserved order parameter
#Time stepping loop
while True:
    VOld.vector()[:] = VNew.vector() #updating the solution of previous time step
    El = assemble(elasticE(eta)*dx)
    Bl = assemble(f1*dx)
    Total_Energy = (El+Bl+assemble((k/2*dot(grad(eta),grad(eta)))*dx))
    TE = [t,El,Bl,Total_Energy-El-Bl,Total_Energy]
    t += dt #updatting dt after each time step
    mu_mid = (1.0-theta)*mu_0 + theta*mu #cranck-nicholson scheme
#variational formulation of the coupled cahn-hilliard equation       
    cahn1 = (eta - eta_0)/dt*V1*dx + M*dot(grad(mu_mid),grad(V1))*dx
    cahn2 = mu*V2*dx -diff(elasticE(eta),(eta))*V2*dx-dfdn*V2*dx -k*dot(grad(eta),grad(V2))*dx
    cahn = cahn1+cahn2
    a = derivative(cahn,VNew,VTrial) # finding the jacobian of the combined non-linear form(cahn)
    problem = CahnHilliard(a,cahn)
    (no_of_iterations,converged) = solver.solve(problem,VNew.vector()) #solving the non-linear equation; no of iterations and whether converged or not are taken as output
    if no_of_iterations <4:  #for adaptive time stepping
        dt = 2*dt
    else:
        dt = dt/2
    file1 << (VNew.split()[0], t) # storing th esolution at each time step
    q = assemble(abs(mu - mu_0)*dx)/len(mesh.coordinates()) #finding the changes in chemical potential wrt to previous time step
    with open('EnergyN_dt_p001.csv',mode='a') as csv_file: # storing the elastic,bulk,grad,toatal energy @each time-step in csv file format
        writer = csv.writer(csv_file)
        writer.writerow(TE)
    if (q < 1e-10): #time-stepping loop terminating condition
        break
print("Loop terminated as the changes in chemical potential reaches to the very small value")
