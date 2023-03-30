close all
clear
clc

% Matrices 
A = [-0.018223  -0.088571   -9.78   0;
    -0.003038   -1.2563     0       1;
    0           0           0       1;
    0.0617      -28.075     0       -4.5937];

B = [0      1.1962;
    0       -0.00120;
    0       0;
    7.84    -4.05];

C = [1  0       0       0;
    0   57.296  0       0;
    0   0       57.296  0;
    0   0       0       57.296;
    0   -57.296 57.296  0];

D = [0  0;
    0   0;
    0   0;
    0   0;
    0   0];

% Obtenir mes valeurs propres
eigA = eig(A)
Ve = ss(A, B, C, D)
% [wnvalid, zetavalid, pvalid] = damp(Ve)
wn = abs(eigA);   %Absolut des poles
zeta = -real(eigA)./wn;
ts = (4)./(zeta.*wn);
Mp = 100.*exp(-pi./tan(acos(zeta)));
wa = wn.*sqrt(1-zeta.^2);
tp = pi./wa; 
FT = tf(Ve)




figure()
hold on
rlocus(FT(1,2))

% Kv = rlocfind(FT(1,2))
Kv = 1.0263;
legend
hold off