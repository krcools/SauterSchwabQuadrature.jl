using CompScienceMeshes
using Test
using SauterSchwabQuadrature

include("verificationintegral.jl")

pI = point(1,5,3)
pII = point(2,5,3)
pIII = point(7,1,0)
pIV = point(5,1,-3)
pV = point(0,0,0)

Sourcechart = simplex(pI,pIII,pII)
Testchart = simplex(pI,pIV,pV)

Accuracy = 12
cv = CommonVertex(Accuracy)

function integrand(x,y)
			return(((x-pI)'*(y-pV))*exp(-im*1*norm(x-y))/(4pi*norm(x-y)))
 end

 function INTEGRAND(û,v̂)
 	n1 = neighborhood(Testchart, û)
 	n2 = neighborhood(Sourcechart, v̂)
 	x = cartesian(n1)
 	y = cartesian(n2)
 	output = integrand(x,y)*jacobian(n1)*jacobian(n2)
  	return(output)
 end

 result = sauterschwab_parameterized(INTEGRAND, cv)-
 			verifintegral2(Sourcechart, Testchart, integrand, Accuracy)

@test norm(result) < 1.e-3
