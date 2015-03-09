
using ecm

begin

	N = BigInt(ARGS[1])

	if isprime(N)
		println(N)
		return
	end

	racine = BigInt(iround(BigFloat(sqrt(N))))

	B = BigInt(iround(BigFloat(exp((1 / sqrt(2)) * sqrt(log(racine) * log(log(racine)))))))

	x = BigInt(1)

	#liste_x = BigInt[]

	for q = 2:B
		if isprime(q)

			x *= q^(iround(log(q, N)))
			#push!(liste_x, x)
		end
	end

	facteur = BigInt(1)
	nb_courbes = 0

	for i = 1:450

		print(nb_courbes)
		print(" ")
		facteur = ecm_f(x, N, 'e')
		nb_courbes += 1

		if facteur != 1 && facteur != N
		
			side = BigInt(N / facteur)
			print(facteur)
			print(" * ")
			println(side)
			return
		end
	end
end