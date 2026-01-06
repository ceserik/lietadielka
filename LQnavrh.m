%3 affine comvination of 3 LQR points

Q = [
    0.1 0 0 0;
    0 0.1 0 0;
    0 0 0.3 0;
    0 0 0  1];
R = 100 ;

rychlosti = [];
vysky = [];
Gainy = [];
for i = 1:3
    FLYmodel
    longK = lqr(podelny_system,Q,R);

    rychlosti(end+1) = plane.v0

    vysky(end+1) = plane.h

    Gainy(end+1,:) = longK

end

% vypocitaj gain scheduling

A = [ones(size(rychlosti))' rychlosti' vysky'];

coeficienty = A\Gainy



