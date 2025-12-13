podelny_system = ss(Ap,Bp,Cp,Dp);
stranovy_system = ss(As,Bs,Cs,Ds);

%% prenos elevatoru na pitch rate
%sisotool(podelny_system(3,1))



%% As = A stranovy %%%%%%%%%%%%%%%%%%%%% STRANOVY SYSTEM
% controls 1.aileron
% controls 2. rudder
% states 1. side slip
% states 2. roll rate
% states 3. yaw rate
% states 4. roll angle
% states 5. yaw angle
% states 6. Accel y




%% legenda
labels = arrayfun(@(i) ...
    sprintf('Speed = %d, Height = %d', legends(i,1), legends(i,2)), ...
    1:size(legends,1), 'UniformOutput', false);





% prenos ailerou na roll rate
%sisotool(stranovy_system(2,1))
figure(2)
bode(stranovy_system(2,1))
hold on
title("Prenos aileron na roll rate")
legend(labels,"Location","best");
grid on
hline = findall(gcf, 'type', 'line');
% Set the linewidth
set(hline, 'LineWidth', 2); % Change 2 to your desired linewidth



% prenos rudder na yaw rate
figure(3)
bode(stranovy_system(3,2))
hold on
title("Prenos rudder na yaw rate")
legend(labels,"Location","best");
grid on
hline = findall(gcf, 'type', 'line');
% Set the linewidth
set(hline, 'LineWidth', 2); % Change 2 to your desired linewidth



%%  Ap = A podelni ############## Podelny system
% controls 1.elevator
% states 1.velocity
% states 2.angle of attack
% states 3.pitch rate
% states 4.accel x
% states 5. accel y


% prenos elevator na pitch rate
figure(4)
bode(podelny_system(3,1))
hold on
title("Prenos rudder na yaw rate")
legend(labels,"Location","best");
grid on

hline = findall(gcf, 'type', 'line');
% Set the linewidth
set(hline, 'LineWidth', 2); % Change 2 to your desired linewidth






