
## Liste des facteurs premiers inférieurs à 30 = 2 * 3 * 5  ##
Facteurs_Premiers_30 = getindex(Int8, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29)

## Liste des facteurs premiers ne divisant pas 30 (+ le facteur 1) ##
Facteurs_speciaux = getindex(Int8, 1, 7, 11, 13, 17, 19, 23, 29)

## -> taux de facteurs testés ~= 26,7 % ~= (10 - 3 + 1) / 30 #

## Factorisation par la méthode naive ##
function factorisation_naive(x)
	
	if(x < 0)
		error("x doit etre positif")
	end

	if isprime(x)
		return x
	end

	# on teste les facteurs premiers inférieurs à 30 #

	for f in Facteurs_Premiers_30

		r = mod(x, f)
		if r == 0 # on a trouvé un facteur
			return f
		end
	end

	## Code à exécuter si le nombre à tester est supérieur à 30^2 ##
		
	arret2::Uint64 = isqrt(x) + 1

	for k = 1:arret2	
		for fs in Facteurs_speciaux
	
			ft::Int64 = 30 * k + fs

			r = mod(x, ft)
			if r == 0 # on a trouvé un facteur
				return ft
			end
		end
	end
end

begin

	N = parseint(ARGS[1])	
	facteur = factorisation_naive(N)
	side = N / facteur
	side = convert(Int64, side)
	println("$facteur * $side")
end


	