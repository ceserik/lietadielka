%3 affine comvination of 3 LQR points
clear
Q = [ 0.1 0 0 0; %speed 
      0 0.7 0 0; %angle of attack 
      0 0 0.05 0; % pitch rate 
      0 0 0 0.05]; % pitch 
R = 200 ;

rychlosti = [];
vysky = [];
Gainy = [];
for i = 1:1
    FLYmodel
    longK = lqr(podelny_system,Q,R);

    rychlosti(end+1) = plane.v0

    vysky(end+1) = plane.h

    Gainy(end+1,:) = longK

end

% vypocitaj gain scheduling

A = [ones(size(rychlosti))' rychlosti' vysky'];

coeficienty = A\Gainy



