clc
close all
clear

%% Variables

% Dynamique Avion
OmegaN = 1;
Zeta = 0.1;
numDynAvion = [(OmegaN.^2)];
denDynAvion = [1 (2.*Zeta.*OmegaN) (OmegaN.^2)];
tfDynAvion = tf(numDynAvion, denDynAvion);

% Compensateur PD
Kp = 1;
Td = [0.0 0.2 0.5 1.0 2.5 5.0];
denComp = [1];

% Actionneur
tau = [0.0 0.01 0.05 0.1 0.2];
numAction = [1];
% denAction = [tau(i) 1];
% tfAction1 = tf(numAction, denAction1);


% Capteur
Tr = [0.05 0.10 0.20 0.25];


% Affichage Liieu des racines
figure('Name','Lieu des racines')
hold on
% Compensateur PD
for i = 1:6
    numComp = [Td(i) 1].*Kp;
    Vc = tf(numComp, denComp);
    sysComp(i) = Vc.*tfDynAvion;
    rlocus(sysComp(i))
    text = convertCharsToStrings(['Td = ', num2str(Td(i))]);
    legendtxt2(i) = [text];
    
end
title('Lieu des racines')
hold off

% 2e graphique
figure('Name','Reponse a un echelon')
hold on
% Compensateur PD
for i = 1:6
    [Y, T] = step(feedback(sysComp(i),1),60);
    plot(T,Y)
end
title('Reponse a un echelon')
xlabel('Time')
ylabel('Amplitude')
legend(legendtxt2)
hold off
grid on
