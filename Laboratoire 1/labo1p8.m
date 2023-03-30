clc
close all
clear

%% Variables
Gcpg = 1;

numGcAvPh = [1 2.9].*4.68;
denGcAvPh = [1 5.4];
GcAvPh = tf(numGcAvPh,denGcAvPh);

numGcPD = numGcAvPh;
denGcPD = [5.4];
GcPD = tf(numGcPD, denGcPD);

numsystem = [4];
densystem = [1 2 0];
system = tf(numsystem, densystem);

%% A
disp('==========A==========')
system1 = Gcpg*system;
system1BF = system1;
p1 = rlocus(system1BF,1);
wn1 = abs(p1);
zeta1 = -real(p1)./wn1;
ts1 = (4)./(zeta1.*wn1);
disp(['Wn1   =  ',num2str(wn1(1))])
disp(['Zeta1 =  ',num2str(zeta1(1))])
disp(['ts1   =  ',num2str(ts1(1))])
%% B
disp('==========B==========')
system2 = GcAvPh*system;
p2 = rlocus(system2,1);
wn2 = abs(p2);
zeta2 = -real(p2)./wn2;
ts2 = (4)./(zeta2.*wn2);
disp(['Zeta2 = ',num2str(zeta2(1))])
disp(['ts2   = ',num2str(ts2(1))])
%% C
disp('==========C==========')
system3 = GcPD*system;
p3 = rlocus(system3,1);
wn3 = abs(p3);
zeta3 = -real(p3)./wn3;
ts3 = (4)./(zeta3.*wn3);
disp(['Zeta3 = ',num2str(zeta3(1))])
disp(['ts3   = ',num2str(ts3(1))])
%% Affichage
figure('Name','Lieu de racines')
hold on
rlocus(system1)
rlocus(system2)
rlocus(system3)
legend('system1','system2','system3')
hold off
