import math
phi = 4.0*math.pi/5.0
sum_cos_col = 0
sum_sin_col = 0

sum_col = 0
for a in range(5):
    sum_col += math.cos(a*phi)
sum_sin = 0
for b in range(5):
    sum_sin+= math.sin(b*phi)

cross_term = 0
for k in range(5):
    cross_term+=math.cos(k*phi)*math.sin(k*phi)
for i in range(5):
    sum_cos_col=math.cos(i*phi)**2

for j in range(5):
    sum_sin_col=math.sin(j*phi)**2

print sum_col
print sum_sin
print cross_term
print sum_cos_col
print sum_sin_col
print 5*((-math.cos(phi)))

print "---"
print "cos   sin"
for a in range(5):
    print math.cos(a*phi), " ", math.sin(a*phi)
