
using ecm

begin

	N = BigInt(ARGS[1])

	#if isprime(N)
		#println(N)
		#return
	#end

	#p = BigInt("1125899906842679")
	p = BigInt("1152921504606847009")
	P = p + 1 + 2 * iround(sqrt(p))
	B = BigInt(iround(BigFloat(exp((1 / sqrt(2)) * sqrt(log(P) * log(log(P)))))))

	x = BigInt(1)

	for q::BigInt = 2:B
		if isprime(q)

			x *= q^(iround(log(q, P)))
		end
	end

	#x = BigInt(1034)
	facteur = BigInt(1)
	nb_courbes = 0

	tic()
	for i = 1:250

		print("$nb_courbes ")
		facteur = ecm_f(x, N, 'e')
		nb_courbes += 1

		if facteur != 1 && facteur != N
		
			toc()
			side::BigInt = N / facteur
			println("$facteur * $side")
			return
		end
	end
end