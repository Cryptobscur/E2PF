
### ECM est pour l'instant configuré sur les courbes de Weierstrass, l'algo sera étendu
### à d'autres courbes lorsque la fiabilité des tests avec Weierstrass sera établie.

type point # Un point en coordonnees projectives #
	X::BigInt
	Y::BigInt
	Z::BigInt
end

point() = point(0, 0, 0)

## Addition de deux points non égaux P et Q sur une courbe de Weierstrass : y^2 = x^3 + ax + b ##
function add_weierstrass(P::point, Q::point, N::BigInt)
	
	R = point()

	## TODO : tester les points particuliers ##

    # on effectue quelques calculs préliminaires pour accélérer le calcul global #
	tmpU = ((Q.Y * P.Z) - (P.Y * Q.Z)) % N
	tmpV = ((Q.X * P.Z) - (P.X * Q.Z)) % N

	tmpUU = tmpU^2 % N # en fait, tmpUU n'est utilisée qu'une fois, mais cela améliore la compréhension
	tmpVV = tmpV^2 % N

	tmpZP_ZQ = (P.Z * Q.Z) % N

	# tmpX contient la valeur R.X / tmpV, elle est utilisée pour le calcul de R.Y #
	tmpX = ((tmpUU * tmpZP_ZQ) + (tmpV * tmpVV)) % N # formule initiale : (tmpU^2 * P.Z * Q.Z) - (((P.X * Q.Z) - (Q.X * P.Z)) * tmpV^2)

	# on calcule les valeurs du point R, somme de P et Q #
	R.Y = ((tmpVV * (tmpU * P.X - tmpV * P.Y) * Q.Z) - (tmpU * tmpX)) % N # formule initiale : (tmpV^2 * ((tmpU * P.X) - (tmpV * P.Y)) * Q.Z) - (tmpU * (R.X / tmpV))
	R.Z = (tmpV * tmpVV * tmpZP_ZQ) % N # formule initiale : tmpV^3 * P.Z * Q.Z
	R.X = (tmpX * tmpV) % N # formule initiale : tmpU^2 * tmpV * P.Z * Q.Z - ((P.X * Q.Z) - (Q.X * P.Z)) * tmpV^3

	return R

end

## Doublement d'un point P (= 2P) sur une courbe de Weierstrass : y^2 = x^3 + ax + b ##
function double_weierstrass(P::point, a::BigInt, N::BigInt)
	
	R = point()

	## TODO : tester les points particuliers ##

    # on effectue quelques calculs préliminaires pour accélérer le calcul global #
	tmpU::BigInt = ((3 * P.X) + (a * P.Z^2)) % N
	tmpV::BigInt = (2 * P.Y * P.Z) % N

	tmpUU::BigInt = tmpU^2 % N # en fait, tmpUU n'est utilisée qu'une fois, mais cela améliore un peu la compréhension
	tmpVV::BigInt = tmpV^2 % N

	tmpVVV::BigInt = (tmpV * tmpVV) % N

	# tmpX contient la valeur R.X / tmpV, elle est utilisée pour le calcul de R.Y #
	tmpX::BigInt = ((tmpUU * P.Z) - (2 * P.X * tmpVV)) % N # formule initiale : (tmpU^2 * P.Z) - (2 * P.X * tmpV^2)

	# on calcule les valeurs du point R, somme de P et Q #
	R.Y = ((tmpU * (P.X * tmpVV - tmpX)) - (tmpVVV * P.Y)) % N # formule initiale : tmpU * (P.X * tmpV^2 - (R.X / tmpV)) - (tmpV^3 * P.Y)
	R.Z = (tmpVVV * P.Z) % N # formule initiale : tmpV^3 * P.Z
	R.X = (tmpX * tmpV) % N # formule initiale : (tmpU^2 * tmpV * P.Z) - (2 * P.X * tmpV^3)

	return R

end

## Algorithme Double & Add pour calculer des expressions du type x * P, avec x un entier et P un point sur une courbe ##
function double_and_add(P::point, x::BigInt, a::BigInt, N::BigInt)
	
	# utile ? #
	if x < 0
		x = -x
		P.Y = -P.Y
	end

	bin_x = digits(x, 2) # on récupère x en base deux
	double_P = point()

	for i = (length(bin_x) - 1):-1:1
		
		double_P = double_weierstrass(P, a, N) # double

		if bin_x[i] == 1
			P = add_weierstrass(P, double_P, N) # add
		end

	end

	return P

end

## à compléter ##
function ecm()

end

