close all
clear all
clc

num = [17.8885];
den = [1 2 0];
Gs = tf(num, den);

[GM, PM, Wp, Wg] = margin(num,den)

disp(['La marge de gain est ', num2str(GM),' à la fréquence ', num2str(Wp), ' rad/s']);
disp(['La marge de phase est ', num2str(PM),' deg à la fréquence ', num2str(Wg), ' rad/s']);
 

% B
num2 = [1 2.5483].*1.5697;
den2 = [1 6.2787];
Ga = tf(num2, den2);
Gs2 = Gs*Ga;

figure('Name','Lieu des racines')
hold on
rlocus(Gs);
rlocus(Gs2,'r');
legend('Original','Avec compenateur')
hold off

figure('Name','Diagramme de Bode')
hold on
bode(Gs);
bode(Gs2,'r');
legend('Original','Avec compenateur')
logspace(-1, 2);
hold off
grid on

% c
figure('Name','Diagramme de Bode de Ga')
hold on
bode(Ga);
legend('Ga')
logspace(-1, 2);
title('Diagramme de Bode du compensateur')
hold off
grid on

