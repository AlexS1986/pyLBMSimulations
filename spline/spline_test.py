import SolidLBM.solver.interpolation.Spline as Spline_Modul
import matplotlib.pyplot as plt

x = list([1,7,12,15,19])
y = list([8,10,7,8,7])

s = Spline_Modul.Spline(x)
s.compute_spline_coefficients(y)



b_test = s.compute_b_from_coefficients_of_derivatives(y)

y1 = s.compute_value_of_spline(x=1.0)
print("hi")

xplt = range(1,19)
yplt = list()
for xx in xplt:
    yy = s.compute_value_of_spline(xx)
    yplt.append(yy)

plt.plot(xplt,yplt)
plt.axis([1, 19, 0, 12])
plt.show()
