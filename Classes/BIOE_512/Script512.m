kon = 5E-6;
koff = .1;
Rtot = .1E-6;
kphos = 1E-6;
kdephos = .01;
Enztot = 10E-6;
GF = 100E-12:1E-8:100E-6;

GFR = (Rtot.*GF)./(GF+(koff/kon));

EnzP = (kphos*Enztot*GFR)./(kdephos+(kphos*GFR));

plot(GF,EnzP);

plot(GF,GFR);