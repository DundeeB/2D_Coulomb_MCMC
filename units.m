% cgs, K=1 (F=q^2/r^2), Kelvin = 1

cm = 1e8;
g = 1;
s = 1;

meter = 100*cm;
Angstrem = 1e-10*meter;
nm = 1e-9*meter;
pm = 1e-12*meter;
r_ion = 282.01*pm;  % https://en.wikipedia.org/wiki/Ionic_radius

Kg = 1e3*g;
statC = cm^(3/2)*g^(1/2)*s^(-1);  % https://en.wikipedia.org/wiki/Statcoulomb

e = 4.80320425e-10*statC;  % statcoulombs F=q^2/r for c g s k=1
Kelvin = 1;

J = Kg*meter^2/s^2;
k_B = 1.38064852e-23*J/Kelvin;