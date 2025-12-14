Q = [
    0.1 0 0 0;
    0 0.1 0 0;
    0 0 0.3 0;
    0 0 0  1];
R = 100 

longK = lqr(podelny_system,Q,R)