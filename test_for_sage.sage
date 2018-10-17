stationary_points = lambda f : solve([gi for gi in f.gradient()], f.variables())
f(x,y,z) = (sqrt(10)*sqrt(x**2+6y**2+3*z**2)/(x + 6y + 3z)**2)**2
g=f(x,y,1)
print stationary_points(g)
