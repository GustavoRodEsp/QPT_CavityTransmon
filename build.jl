module build
push!(LOAD_PATH, pwd())
using QuantumOptics

function HamiltonianCT(Nmaxa,Nmaxb,Delta,xi,alpha,k,xip)
   Fa = FockBasis(Nmaxa)
   Fb = FockBasis(Nmaxb)
   a = destroy(Fa)⊗one(Fb)
   ad= dagger(a)
   b = one(Fb)⊗destroy(Fb)
   bd= dagger(b)
   #H=Delta*ad*a + Delta*bd*b
   #H=Delta*ad*a-xi*ad*a*bd*b-(alpha/2)*bd*bd*b*b-(k/2)*ad*ad*a*a-(xip/2)bd*b*ad*ad*a*a
   H=Delta*ad*a-0*(xi*ad*a*bd*b-(alpha/2)*bd*bd*b*b-(k/2)*ad*ad*a*a-(xip/2)bd*b*ad*ad*a*a)
   return H
end


function wignerF(rho,L,N,h,name1,name2)
rhoa=ptrace(rho,1)
  rhob=ptrace(rho,2)
  open(name1,"w") do io
  open(name2,"w") do io2
  dx = 2L/N
  global negcav=0
  global negtrans=0
  global volcav=0
  global voltrans=0
  for x in -N:dx:N
	for p in -N:dx:N
	   wigcav = wigner(rhoa,x,p)
	   wigtrans = wigner(rhob,x,p)
	   global negcav +=dx*dx*(abs(wigcav)) 
	   global negtrans += dx*dx*(abs(wigtrans))
	   global volcav += dx*dx*wigcav
	   global voltrans += dx*dx*wigtrans
	   println(io , x/sqrt(h)," ",p/sqrt(h)," ", wigcav)
	   println(io2 , x/sqrt(h)," ",p/sqrt(h)," ", wigtrans)
	end
  end
  end #De Open2
  end #De Open1
  return [volcav, negcav-1,voltrans, negtrans-1]
end





end
