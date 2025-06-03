module main_stationary
push!(LOAD_PATH, pwd())
using QuantumOptics
import build

#Hamiltonian parameters
L=12
N=200
hh=2
Nmaxa=6
Nmaxb=6
Delta=1
xi=0.000000001
alpha=0.000000001
k=0.000000001
xip=0.000000001
n=1
name1="output/wigrhoa.dat"
name2="output/wigrhob.dat"

h=0.1
deltalist=[i*h for i in -50:50]
open("output/spectrum.dat","w") do io
for i in deltalist
  Ham = build.HamiltonianCT(Nmaxa,Nmaxb,i,xi,alpha,k,xip)
  #Ham = build.HamiltonianCT(Nmaxa,Nmaxb,i,xi,alpha,k,xip)
  eg = eigenstates(dense(Ham))
  println(io , i," ",eg[1][1])
end
end

Ham = build.HamiltonianCT(Nmaxa,Nmaxb,Delta,xi,alpha,k,xip)
esystem= eigenstates(dense(Ham))
psi=esystem[2][n]
rho = dm(psi)
w = build.wignerF(rho,L,N,hh,name1,name2)
println(w)

end #Este end es para cerrar el modulo