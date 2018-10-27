M.<x,y> = PolynomialRing(QQ);
Mx.<X> = PolynomialRing(QQ);

def get_elem (n, d, a, b, L):
	nx = (n/d).derivative(y,b)(y=0).numerator()/factorial(b);
	print "nx: ", nx
	dx = (n/d).derivative(y,b)(y=0).denominator()
	print "dx: ", dx
	nX = Mx(0)
	dX = Mx(0)

	for i in range(nx.degree(x)+1):
		nX += nx.coefficient({x:i, y:0}).constant_coefficient() * X^i
	for i in range(dx.degree(x)+1):
		dX += dx.coefficient({x:i, y:0}).constant_coefficient() * X^i
	
	print "nX: ", nX
	print "dX: ", dX

	res = (nX * dX.inverse_mod(X^(a+L+1))) % X^(a+L+1)
	return res
