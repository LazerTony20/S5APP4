clc
close all
clear
clc

% Valeurs des simulations
t1 = [0:1:500];
t2 = [0:0.1:10];
u1 = ones(size(t1));
u2 = ones(size(t2));
%

A = [-0.018223 -0.088571 -9.78 0;
    -0.003038 -1.2563 0 1;
    0 0 0 1;
    0.0617 -28.075 0 -4.5937];
B = [0 1.1962;
    0 -0.00120;
    0 0;
    7.84 -4.05];
C = [1 0 0 0;
    0 57.296 0 0;
    0 0 57.296 0;
    0 0 0 57.296;
    0 -57.296 57.296 0];
D = [0 0;
    0 0;
    0 0;
    0 0;
    0 0];
ss0 = ss(A,B,C,D);
disp(['Le systeme de variable d''état obtenu avec ss: ']), ss0
%% A 
%analyse des caractéristiques dynamiques de l'avion à partir des valeurs
%propres du système et de la réponse temporelle de la vitesse v du système
%soumis à un échelon sur a_p : donner le temps du premier pic, le
%dépassement maximum, le temps de stabilisation, le facteur
%d'amortissement, la période des oscillations amorties et naturelles et
%vérifier ces résultats à partir de la réponse temporelle
p=eig(A);
wn = abs(p);
zeta = real(p)./wn;
ts = -4./(zeta.*wn);
disp(['ts initial = ' num2str(ts(3))])
wa = wn.*sqrt(1-zeta.^2);
tp = pi./wa;
phi = acos(zeta);
mp = 100 * exp(-pi./tan(phi))

FTBO = tf(ss0)

figure('Name','Mandat A')
hold on
plot(t1,lsim(FTBO(1,2),u1,t1))
yline(mp(3), 'Mp')
xline(tp(3))
xline(ts(3))
legend()
grid on
hold off

%% B
%identification de la fonction de transfert à phase non-minimale à partir
%des pôles et zéros

figure('Name','Determination des fonctions a phase non-minimale')
hold on
subplot(5,2,1)
pzmap(FTBO(1,1))
subplot(5,2,2)
pzmap(FTBO(2,1))
subplot(5,2,3)
pzmap(FTBO(3,1))
subplot(5,2,4)
pzmap(FTBO(4,1))
subplot(5,2,5)
pzmap(FTBO(5,1))
subplot(5,2,6)
pzmap(FTBO(1,2))
subplot(5,2,7)
pzmap(FTBO(2,2))
subplot(5,2,8)
pzmap(FTBO(3,2))
subplot(5,2,9)
pzmap(FTBO(4,2))
subplot(5,2,10)
pzmap(FTBO(5,2))
hold off
%% C
% présentation du lieu des racines fait à la main entre a_p et v (avec étapes)
%et validation sur MATLAB
%[num,den] = tfdata(FTBO(1,2));
%[P,Z]=pzmap(FTBO(1,2));
figure('Name','Comparaison Rlocus avec à la main')
rlocus(FTBO(1,2))
%% D
%explication de l'effet de la rétroaction Kv sur la stabilité, le temps de
%stabilisation et le dépassement maximum à partir du lieu des racines
% Lorsqu'on change le gain de la boucle de rétroaction,
% les poles se déplacent sur les branches du lieu de
% racine. Avec de nouveaux poles, nous avons des nouvelles
% valeurs de zeta, temps de stabilisation et pic maximal.
K0 = [0:0.01:4];
R0 = rlocus(FTBO(1,2),K0);
W0=abs(R0);
zeta0 = real(R0)./W0;
ts0 = -4./(zeta0.*W0);
wa0 = W0.*sqrt(1-zeta0.^2);
tp0 = pi./wa0;
phi0 = acos(zeta0);
mp0 = 100 * exp(-pi./tan(phi0));
figure('Name','Optimisation du Kv')
hold on
subplot(2,2,1)
plot(K0,ts0)
title('Impact de Kv sur le temps de stabilisation')
subplot(2,2,2)
plot(K0,wa0)
title('Impact de Kv sur la stabilite')
subplot(2,2,3)
plot(K0,phi0)
title('Impact de Kv sur la phase')
subplot(2,2,4)
plot(K0,mp0)
title('Impact de Kv sur le pic maximal')
hold off
%% E
% conception de la boucle interne à partir du lieu des racines: 
%(1) choix de kv
%(2) calcul du nouveau modele n1/d1 et du modele A1, B1, C1, D1 incluant
%la boucle interne
%On voit que dans le mode phigoide, il faut un gain Kv de 1.0263 pour minimiser
%w*z bien que les autres poles seront également légerement déplacé, cet
%impact sera beaucoup moindre que l'impact sur le mode phigoide

[NumFTo4, DenFTo4] = tfdata(FTBO(1,2),'v');
NumDiffDen = conv([0 3.*NumFTo4(2) 2.*NumFTo4(3) NumFTo4(4)],DenFTo4);
NumDenDiff = conv(NumFTo4,[4.*DenFTo4(1) 3.*DenFTo4(2) 2.*DenFTo4(3) DenFTo4(4)]);
NpDNDp = NumDiffDen - NumDenDiff;
FTBOJunc = roots(NpDNDp);
sjonction = FTBOJunc(5);
Kv = -(DenFTo4(1).*sjonction.^4 + DenFTo4(2).*sjonction.^3 + DenFTo4(3).*sjonction.^2 + DenFTo4(4).*sjonction + DenFTo4(5))./(NumFTo4(2).*sjonction.^3 + NumFTo4(3).*sjonction.^2 + NumFTo4(4).*sjonction + NumFTo4(5))

disp(['Nouvelle fonction de transfert avec la retroaction Kv=', num2str(Kv)] )
% Kv=1.0263;  %Valeur obtenue avec rlocfind sur mon lieu de racine de ftbo, entree a bras pour simplifier la vie
FTBF = feedback(FTBO(1,2),Kv)
sim1=lsim(FTBO(1,2),u1,t1);
sim2=lsim(FTBF,u2,t2);
figure('Name','Comparaison avant-apres feedback')
hold on
subplot(2,1,1)
plot(t1,sim1)
title('Stabilisation originale')
grid on
subplot(2,1,2)
plot(t2,sim2)
title('Stabilisation avec retroaction Kv')
grid on
disp('Creation des nouvelles matrice A1,B1,C1,D1')

ss1 = feedback(ss0,[0 Kv; 0 0; 0 0; 0 0; 0 0]')
A1=ss1.A
B1=ss1.B
C1=ss1.C
D1=ss1.D
FTBO2 = tf(ss1)
disp('Calcul du nouveau Ts:')
%[numf,denf]=tfdata(FTBF,'v');
%[rf,pf,kf] = residue(numf,denf)
pf=eig(A1)
wnF = abs(pf)
zetaF = real(pf)./wnF;
tsF = -4./(zetaF.*wnF)
waF = wn.*sqrt(1-zetaF.^2);
tpF = pi./waF;
phiF = acos(zetaF);
mpF = 100 * exp(-pi./tan(phiF));
%% F
% vérification des marges avec diagrammes de Bode et Kv:commentaire sur le
%sens des marges, utilité
figure('Name','Bode Plot et marges')
subplot(2,1,1)
margin(FTBO(1,2))
grid on
subplot(2,1,2)
margin(FTBF)
grid on
%[mag,pha]=bode(FTBF)
hold off
%% G
%réduction de la FT entre a_p et v, présentation de la nouvelle FT réduite, fermeture de la boucle à
%la main sur cette FT, explication de l’effet de Kv sur les  paramètres standards
%(K,zeta,wn,tau) et la reponse
[R1,P1,K1]=residue(FTBO(1,2).numerator{1},FTBO(1,2).denominator{1});
Cdom=abs(R1)./abs(real(P1));
[NUMR,DENR]=residue(R1(3:4),P1(3:4),K1);
TFR = tf(NUMR,DENR);
TFRnocomp = TFR;
TFR = (dcgain(FTBO(1,2))./dcgain(TFR))*TFR;
[NUMR,DENR] = tfdata(TFR,'v');
% [ZR, PR, KR] = tf2zp(NUMR,DENR)
figure('Name','Root locus de la fonction reduite')
rlocus(TFR)
title('Root locus de la fonction reduite compensee')
legend('Ordre 2 Reduite')
%Calcul à la main effectué. On regarde maintenant l'impact de Kv avec la
%valeur trouvée en E)
% Kv = 1.298549011;     % Valeur Trouvee a la main comme reference
% On utilise une partie du travail effectue a la main
% NUMRDIFF = [0 NUMR(1)];     %Derivee de NUMR
% DENRDIFF = [2 DENR(2)];     %Derivee de DENR
% Sseparation = roots(conv(NUMR,DENRDIFF) - conv(NUMRDIFF,DENR))
% KvFB = - (DENR(1).*Sseparation(1).^2 + DENR(2).*Sseparation(1) + DENR(3))./(NUMR(1).*Sseparation(1) + NUMR(2))
% KvFBroots = roots([1.340964 -1.44243704 -0.04699232])

KvFBcoef1 = (NUMR(end-1)./2).^2;
KvFBcoef2 = (NUMR(end-1).*DENR(2) - NUMR(end));
KvFBcoef3 = ((DENR(2)./2).^2 - DENR(3));
KvFBroots = roots([KvFBcoef1 KvFBcoef2 KvFBcoef3]);
KvFB = KvFBroots(1);
NUMR2 = KvFB.*NUMR;
DENR2 = DENR + [0 KvFB.*NUMR(end-1) KvFB.*NUMR(end)];
TFRfb = tf(NUMR2,DENR2)
% pr2 = rlocus(TFRfb,1)
% wnr2 = abs(pr2)
% zetar2 = -real(pr2)./wnr2
[WnTFRfb, ZetarTFRfb, pTFRfb] = damp(TFRfb)
% zerosR2 = roots(NUMR2)
tsr2 = (4)./(ZetarTFRfb.*WnTFRfb)

figure('Name','Reponse a un echelon de la nouvelle fonction de transfert reduite avec boucle fermee kv')
% plot(step(TFRfb,10))
hold on
subplot(2,1,1)
hold on
plot(t1,sim1)
plot(t1,lsim(TFR,u1,t1))
grid on
legend('Ordre 4','Ordre 2 Reduit')
title('Stabilisation Ordre 4 vs Ordre 2 Reduit')
subplot(2,1,2)
hold on
plot(t2,sim2)
plot(t2,lsim(TFRfb,u2,t2))
title('Stabilisation retroaction Kv')
legend(['Ordre 4 Kv=' num2str(Kv)],['Ordre 2 Kv=' num2str(KvFB)])
grid on
hold off
%Sur ce nouveau lieu de racine, nous remarquons que 
%% H
% Calcul du Kv avec le système réduit et comparaison avec la valeur trouvée ci-dessus en (e)
%% I
% comparaison des lieux des racines original et réduit de la FT entre a_p
% et v sur le meme graphique
figure('Name','Comparaison des root locus')
hold on
% subplot(2,1,1)
rlocus(FTBO(1,2),TFRnocomp,TFR)
legend('FTBO Ordre 4','FTBO Ordre 2 reduit','FTBO Ordre 2 Reduit comp DC')
% subplot(2,1,2)
% pzmap(TFRfb)
% legend('FTBF Reduit')
xlim([-3 0.5])
% ylim([-1 1])
 ylim([-8 8])
hold off

% figure()
% hold on
% rlocus(FTBO(1,2),TFR)
% hold off
% legend()
%% J
% Lieu des racines entre gama_d et gama incluant la boucle interne (i.e.
% n1/d1), effet de Kp sur la reponse
figure('Name','Lieu des racines entre gama_d et gama incluant la boucle interne')
hold on
subplot(2,1,1)
rlocus(FTBO(5,1))
subplot(2,1,2)
rlocus(FTBO2(5,1))

%% K
% design de la boucle externe: avec Bode,
% (1)calculer Kp pour instabilite
% (2)calculer Kp pour marges de 6db et 30deg
% (3)mesurer l'erreur en regime permanent par simulation sur MATLAB
% [numKp, denKp] = pade(FTBO2(5,1),4)
% sgrid(0,[]); % specify the desired location on the real axis
% KpInst = rlocfind(FTBO2(5,1))
[NumKpInst, DenKpInst] = tfdata(FTBO2(5,1),'v')

KpInst = DenKpInst(5)./NumKpInst(5)
% KpInst = 0.278;

FTBF2 = feedback(FTBO2(5,1),KpInst);

sim3=lsim(FTBO2(5,1),u2,t2);
sim4=lsim(FTBF2,u2,t2);

figure('Name','Lieu des racines avec Kp pour Instabilite')
hold on
subplot(2,1,1)
plot(t2,sim3)
title('FTBO')
grid on
subplot(2,1,2)
plot(t2,sim4)
title('FTBF instable with Kp')
grid on



%% L
% calcul de l'erreur en régime permanent à partir du gain statique tel que lu sur le diagramme de Bode,
%comparaison avec (3) de l’item ci-dessus
%% M
%calcul du nouveau modèle n2/d2 et le modele A2, B2, C2, D2
%% N
% comparaison et discussion des réponses temporelles avec compensateurs P, PD, PI et PID