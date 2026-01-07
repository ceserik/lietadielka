% Given data
Speed  = [73; 186; 256];      % speed
Height = [143; 12633; 224];   % height
Gain   = [0.68; 0.0287; 0.0309];           % P gains

% Build regression matrix
A = [ones(size(Speed)) Speed Height];

% Solve for coefficients (exact fit for 3 points)
coeff = A \ Gain;

% Extract coefficients
a0 = coeff(1);
a1 = coeff(2);
a2 = coeff(3);

% Display result
fprintf('K = %.6f + %.6f*Speed + %.6f*Height\n', a0, a1, a2);
