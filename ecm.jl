
module ecm

export ecm_f, double_and_add

type point # Un point en coordonnees projectives #
	X::BigInt
	Y::BigInt
	Z::BigInt
end

### WEIRSTRASS ###

## Addition de deux points non égaux P et Q sur une courbe de Weierstrass : y^2 = x^3 + ax + b ##
function add_weierstrass(P::point, Q::point, a::BigInt, N::BigInt)
	
	R = point(0, 1, 0)

	## Si un des points est le point à l'infini (i.e l'élément neutre)
	## alors la somme aura pour résultat l'autre point 
	
	if P.Z == 0
		return Q
	elseif Q.Z == 0
		return P
	end

    # on effectue quelques calculs préliminaires pour accélérer le calcul global #

	tmpV::BigInt = mod((Q.X * P.Z) - (P.X * Q.Z), N)
	tmpU::BigInt = mod((Q.Y * P.Z) - (P.Y * Q.Z), N)

	if tmpV == 0 # alors P.X = Q.X
		if tmpU == 0 # alors P.Y = Q.Y donc P = Q
			return double_weierstrass(P, a, N)
		else # Q = -P donc on additionne P et -P
			return point(0, 1, 0)
		end
	end
	
	tmpUU::BigInt = mod(tmpU^2, N)
	tmpVV::BigInt = mod(tmpV^2, N)

	tmpZP_ZQ::BigInt = mod(P.Z * Q.Z, N)

	# tmpX contient la valeur R.X / tmpV, elle est utilisée pour le calcul de R.Y #

	tmpX::BigInt = mod((tmpUU * tmpZP_ZQ) - ((P.X * Q.Z + Q.X * P.Z) * tmpVV), N) # formule initiale : (tmpU^2 * P.Z * Q.Z) - (((P.X * Q.Z) - (Q.X * P.Z)) * tmpV^2)

	# on calcule les valeurs du point R, somme de P et Q #

	R.Y = mod((tmpVV * (tmpU * P.X - tmpV * P.Y) * Q.Z) - (tmpU * tmpX), N) # formule initiale : (tmpV^2 * ((tmpU * P.X) - (tmpV * P.Y)) * Q.Z) - (tmpU * (R.X / tmpV))
	R.Z = mod(tmpV * tmpVV * tmpZP_ZQ, N) # formule initiale : tmpV^3 * P.Z * Q.Z
	R.X = mod(tmpX * tmpV, N) # formule initiale : tmpU^2 * tmpV * P.Z * Q.Z - ((P.X * Q.Z) - (Q.X * P.Z)) * tmpV^3

	return R

end

## Doublement d'un point P (= 2P) sur une courbe de Weierstrass : y^2 = x^3 + ax + b ##
function double_weierstrass(P::point, a::BigInt, N::BigInt)
	
	R = point(0, 1, 0)

	# on effectue quelques calculs préliminaires pour accélérer le calcul global #

	## Si P est un point 2-torsion alors 2P = 0 (point à l'infini), ou encore y = 0 (ne pas confondre avec P.Y)
	## Si P est le point à l'infini (l'élément neutre), le multiplier par deux donnera 
	## également le point à l'infini, on s'évite donc les tests exécutés dans l'addition.

	tmpV::BigInt = mod(2 * P.Y * P.Z, N)

	if tmpV == 0 # alors on a y = 0 (P.Y * P.Z = (y/z) * z = y)
		return point(0, 1, 0)
	end
	
	tmpU::BigInt = mod(3 * (P.X)^2 + a * (P.Z)^2, N)

	tmpUU::BigInt = mod(tmpU^2, N)
	tmpVV::BigInt = mod(tmpV^2, N)

	tmpVVV::BigInt = mod(tmpV * tmpVV, N)

	# tmpX contient la valeur R.X / tmpV, elle est utilisée pour le calcul de R.Y #

	tmpX::BigInt = mod((tmpUU * P.Z) - (2 * P.X * tmpVV), N) # formule initiale : (tmpU^2 * P.Z) - (2 * P.X * tmpV^2)

	# on calcule les valeurs du point R, somme de P et P (= 2P) #

	R.Y = mod((tmpU * (P.X * tmpVV - tmpX)) - (tmpVVV * P.Y), N) # formule initiale : tmpU * (P.X * tmpV^2 - (R.X / tmpV)) - (tmpV^3 * P.Y)
	R.Z = mod(tmpVVV * P.Z, N) # formule initiale : tmpV^3 * P.Z
	R.X = mod(tmpX * tmpV, N) # formule initiale : (tmpU^2 * tmpV * P.Z) - (2 * P.X * tmpV^3)

	return R

end

### EDWARDS ###

## Addition de deux points non égaux P et Q sur une courbe d'Edwards : x^2 + y^2 = 1 + d * x^2 * y^2 ##
function add_edwards(P::point, Q::point, d::BigInt, N::BigInt)

	R = point(0, 1, 1)

	## le calcul trouve sa source ici : http://www.hyperelliptic.org/EFD/g1p/auto-edwards-projective.html#addition-add-2007-bl ##

	A::BigInt = mod(P.Z * Q.Z, N)
	B::BigInt = mod(A^2, N)

	C::BigInt = mod(P.X * Q.X, N)
	D::BigInt = mod(P.Y * Q.Y, N)

	E::BigInt = mod(d * C * D, N)

	F::BigInt = mod(B - E, N)
	G::BigInt = mod(B + E, N)

	R.X = mod(A * F * ((P.X + P.Y) * (Q.X + Q.Y) - C - D), N)		
	R.Y = mod(A * G * (D - C), N)
	R.Z = mod(F * G, N) # formule initiale : c * F * G

	return R
end

## Doublement d'un point P (= 2P) sur une courbe d'Edwards : x^2 + y^2 = 1 + d * x^2 * y^2 ##
function double_edwards(P::point, d::BigInt, N::BigInt)

	# l'addition est unifiée sur les courbes d'Edwards, i.e pas de cas particulier pour deux points égaux #

	return add_edwards(P, P, d, N);
end

## Algorithme Double & Add pour calculer des expressions du type x * P, avec x un entier et P un point sur une courbe ##
## 'vs' désigne la valeur dite "statique" dans la courbe elliptique employée (une constante en gros) ##
function double_and_add(P::point, x::BigInt, vs::BigInt, N::BigInt, modele::Char) 

	bin_x = digits(x, 2) # on récupère x en base deux
	
	## julia n'a pas encore de switch (officiel du moins) pour l'instant... ##

	if modele == 'w'

		R = point(0, 1, 0)
		double_P = P

		for i = 1:length(bin_x)

			if bin_x[i] == 1
				R = add_weierstrass(R, double_P, vs, N) # add
			end

			double_P = double_weierstrass(double_P, vs, N) # double ##### POURQUOI P EN PARAMETRE ?!!!!! #
		end

	elseif modele == 'e'

		R = point(0, 1, 1)
		double_P = P

		for i = 1:length(bin_x)

			if bin_x[i] == 1
				R = add_edwards(R, double_P, vs, N) # add
			end

			double_P = double_edwards(double_P, vs, N) # double 
		end
	end

	return R
end

## à finaliser ##
function ecm_f(x::BigInt, N::BigInt, modele::Char)

	## julia n'a pas encore de switch (officiel du moins) pour l'instant... ##

	if modele == 'w' # Weirstrass

		## On initialise une courbe et un point sur cette courbe.
		## On fixe la coordonnée Z à 1 afin de se retrouver sous la forme y^2 = x^3 + ax + b
		## De cette façon, b est déterminé par le calcul de y^2 - x^3 - ax
		## et on obtient un point valable (au final, on ne se sert pas de b, donc à ignorer)

		# surement possible de faire qqc de plus "propre" #
		a::BigInt = iround(BigFloat(rand() * (N - 1)))

		## si jamais on fixe b à 1 (bien que cela ne change rien du point de vue du code) on peut alors
		## prendre d'office le point (0:1:1) et on s'évite les deux lignes de code ci-dessous, mais est-ce un choix judicieux ?

		Xw::BigInt = iround(BigFloat(rand() * (N - 1)))
		Yw::BigInt = iround(BigFloat(rand() * (N - 1)))

		Pw = point(0, 1, 1) # 0 et 1 ou Xw et Yw ?

		Rw::point = double_and_add(Pw, x, a, N, 'w')

		println(Rw)
		val::BigInt = BigInt(gcd(Rw.Z, N))
		println("PGCD = $val")
		return val

	elseif modele == 'e' # Edwards

		## On initialise une courbe et un point sur cette courbe.
		## On fixe la coordonnée Z à 1 afin de se retrouver sous la forme x^2 + y^2 = 1 + d * x^2 * y^ 2
		## De cette façon, d est déterminé par le calcul de (x^2 + y^2 - 1) / (x^2 * y^2)
		## et on obtient un point valable (l'inversion est ici inévitable ?)

		Xe::BigInt = iround(BigFloat(rand() * (N - 1)))
		Ye::BigInt = iround(BigFloat(rand() * (N - 1)))

		Xe2 = mod(Xe^2, N)
		Ye2 = mod(Ye^2, N)

		tmp = Xe2 * Ye2 # = x^2 * y^2

		## Julia est plutot "chiant" au niveau de gmp car il transforme les valeurs de retour de certaines fonctions en erreurs/exceptions
		## Il faut donc aller voir le code julia associé pour déterminer le type d'erreur et se coltiner sur le dos des try/catch...

		try tmp = invmod(tmp, N) # on tente d'inverser x^2 * y^2
		catch test
			if isa(test, ErrorException) # si l'inversion n'est pas possible, on a trouvé un facteur de N
				println("inversion impossible !")
				return BigInt(gcd(tmp, N))
			end
		end

		d::BigInt = mod((Xe2 + Ye2 - 1) * tmp, N) # d = (x^2 + y^2 - 1) / (x^2 * y^2)

		Pe = point(Xe, Ye, 1)

		Re::point = double_and_add(Pe, x, d, N, 'e')

		println(Re)
		val2::BigInt = BigInt(gcd(Re.X, N))
		println("PGCD = $val2")
		return val2
	end
end

end