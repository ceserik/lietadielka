function x_trim = trimConfig(plane)
v0 = plane.v0;
m = plane.m;
M = plane.M;
rho = plane.rho;
g = 9.81;

Ixx = plane.Jxx;    % moment setrvacnosti pro osu x
Izz = plane.Jyy;    % moment setrvacnosti pro osu y
Iyy = plane.Jzz;    % moment setrvacnosti pro osu z
Ixz = plane.Jxy;    % moment setrvacnosti pro osu z
cfgmatfile = 'MyAircraft';

% stredni aerodynicka tetiva
c = 2.15; % m
% rozpeti
b = 9.17; % m
% plocha kridla
S = 18.8; % m^2


%% Aerodynamicke charakteristiky
% vztlakovy koeficient
    % pocatecni hodnota vztlakoveho soucinitele
    x_cy0       = [0,0.4,0.6,0.7,0.75,0.8,0.85,0.9]; 
    y_cy0       = [0.09,0.09,0.09,0.105,0.13,0.12,0.07,0.02];
    CL0 = interp1(x_cy0,y_cy0,M);
    
    % derivace soucinitele vztlaku podle uhlu nabehu
    if M <= 0.3
        CLa = 4.4418;
    else
        CLa = - 62.716*M^5 + 122.09*M^4 - 85.797*M^3 + 28.153*M^2 ...
                   - 3.8433*M + 4.541;
    end;

    % derivace soucinitele vztlaku podle vychylky vyskovky
    CLde = 0.305;

    % derivace soucinitele vztlaku podle bezrozmerne zmeny uhlu nabehu
    if M<0.3
        Mp = 0.3;
        CLalphadot = -49.167 + 558.542*Mp - 2284.051*Mp^2 ... 
        + 4412.037*Mp^3 - 4067.130*Mp^4 + 1435.185*Mp^5;
    elseif M>0.75
        Mp = 0.75;
        CLalphadot = -49.167 + 558.542*Mp - 2284.051*Mp^2 ... 
        + 4412.037*Mp^3 - 4067.130*Mp^4 + 1435.185*Mp^5;
    else
        CLalphadot = -49.167 + 558.542*M - 2284.051*M^2 ... 
        + 4412.037*M^3 - 4067.130*M^4 + 1435.185*M^5;
    end

    % derivace soucinitele vztlaku podle bezrozmerne rychlosti kloneni
    CLq = 5.7414 - 9.3115*M + 20.8*M^2 - 13.465*M^3;


% odporovy soucinitel
    % pocatecni hodnota odporoveho soucinitele
    CD0 = 0.08;
    
    % derivace soucinitele odporu podle uhlu nabehu
    if M<0.3
        Mp = 0.3;
        CDa = -1*(10.081 - 23.378*Mp + 13.048*Mp^2);
    else
        CDa = -1*(10.081 - 23.378*M + 13.048*M^2);
    end

    % derivace soucinitele odporu podle vychylky vyskovky
    if M<0.3
        Mp = 0.3;
        CDde =( 0.1857*Mp^2 - 0.1897*Mp + 0.0649);
    elseif M>0.7
        Mp = 0.7;
        CDde =( 0.1857*Mp^2 - 0.1897*Mp + 0.0649);
    else
        CDde =( 0.1857*M^2 - 0.1897*M + 0.0649);
    end


% bocny soucinitel
    cyPom = m*9.81/(0.5*rho*v0^2*S);
    
    % derivace soucinitele bocne sily podle uhlu vyboceni
    CYbeta = (-0.946 - 0.0333*cyPom + 0.0032*M - 0.16*M^2);
    
    % derivace soucinitele bocne sily podle vychylky kridelek
    CYda = -0.02;

    % derivace soucinitele bocne sily podle vychylky smerovky
    CYdr = (-0.1845 + 0.0325*cyPom -0.0093*M -0.0228*M^2 ...
        -0.012*M^3);

    %  derivace soucinitele bocne sily podle kloneni
    CYp = 0;
    
    %  derivace soucinitele bocne sily podle zatacive rychlosti
    CYr = 0;

    
% soucinitel klopiveho momentu
    y_mz0 = [0.025,0.025,0.03,0.025,0.02,0.05,-0.02,-0.03];
    % pocatecni hodnota klopiveho soucinitele
    Cm0 = interp1(x_cy0,y_mz0,M);

    % derivace soucinitele klopiveho momentu podle uhlu nabehu
    if M <= 0.3
        Cma = -0.4899;
    else
        Cma = 5.9964*M^4 - 10.271*M^3 + 6.7921*M^2 - 2.0056*M -0.2708;
    end    

    % derivace soucinitele klopiveho momentu podle vychylky vyskovky
    Cmde = -0.88;

    % derivace soucinitele klopiveho momentu podle zmeny uhlu nabehu
    if M <= 0.3
        Cmalphadot = -1.38271969;
    else
        Cmalphadot = 502.427*M^5 - 1303.11*M^4 + 1330.2481*M^3 ...
            - 668.0225*M^2 + 164.08*M - 17.0671;
    end    
    
    % derivace soucinitele klopiveho momentu podle klopeni
    if M <= 0.3
        Cmq = -4.9934668;
    else
        Cmq = 2.5926*M^3 - 8.7063*M^2 + 4.8590*M - 5.7376;
    end    

    
% soucinitel kloniveho momentu
    % derivace soucinitele klopiveho momentu podle uhlu vyboceni
    if M<0.45
        Clbeta = -0.056 - 0.0576*cyPom - 0.0088*M;
    else
        Clbeta = -0.06 - 0.0576*cyPom - 0.0285*(M - 0.45);
    end
    
    % derivace soucinitele klopiveho momentu podle vychylky kridelek
    if M<=0.6
        Clda = (-0.156 + 0.04*(M -0.2));
    else
        Clda = -1.1039 + 4.6696*M - 7.5941*M^2 ...
            + 4.1483*M^3;
    end
    
    % derivace soucinitele klopiveho momentu podle vychylky smerovky
    Cldr = (-0.03 + 0.008*cyPom)*1;
    
    % derivace soucinitele klopiveho momentu podle kloneni
    Clp = -0.5365 - 1.5267*M + 10.781*M^2 ...
        -20.666*M^3 + 12.4242*M^4;
    
    % derivace soucinitele klopiveho momentu podle zatacive rychlosti
    Clr = (-0.075 - 0.215*cyPom -0.016*M);

% soucinitel zataciveho momentu
    % derivace soucinitele zataciveho momentu podle uhlu vyboceni
    Cnbeta = (-0.173 + 0.01667*cyPom + 0.0092*M -0.1*M^2 +...
        0.0208*M^3);
    
    % derivace soucinitele klopiveho momentu podle vychylky kridelek
    Cnda = -0.0015;

    % derivace soucinitele klopiveho momentu podle vychylky smerovky
    Cndr = (0.015*cyPom - 0.095 -0.0056*M -0.0089*M^2 ...
        - 0.0068*M^3);
    
    % derivace soucinitele klopiveho momentu podle kloneni
    if M<0.7
        Cnp = -0.07435 + 0.045*cyPom - 0.0092*M;
    else
        Cnp = -0.08079 + 0.044*cyPom - 0.0333*(M-0.7);
    end

    % derivace soucinitele klopiveho momentu podle zatacive rychlosti
    Cnr = -0.2291 - 0.0065*cyPom -0.0236*cyPom^2 + 0.0069*cyPom^3 ...
        - 0.0219*M + 0.0097*M^2 - 0.0447*M^3;


%% trimovani letounu

CG = 0.245 ; 
rAC = plane.x_T;
Jgross = [Ixx Iyy Izz Ixz];
qb = 1/2*rho*v0^2;
disp('-----------------------------------')
disp(' Trim of aircraft?')
%vyber = input(sprintf(' [0 - no] [1 - yes]:'));
vyber =1;
x_trim = [ 0 0 0 0];
x_T0 = 0.2450;
x_T = plane.x_T;
if vyber == 1
    if M<0.4
        CD0p = 0.0223;
    elseif M>0.92 
        CD0p = 0.1032;
    else
        CD0p = 7.861*M^4 - 18.378*M^3 + 15.832*M^2 - 5.9269*M + 0.8349;
    end
    qb = 1/2*rho*v0^2;
    fnc = @(x)((Cma*x(1) + Cmde*x(2) + Cm0 + (x_T - x_T0)*(CL0 + CLa*x(1) + CLde*x(2)))^2 +...
        (CLa*x(1) + CLde*x(2) - m*g/(qb*S)*cos(x(3)) + CL0)^2 +...
        (CD0p + 0.00103*(x(1)*180/pi)^2 - 0.00308*x(1)*180/pi +... 
        6.2863E-07*(x(2)*180/pi)^3 + 4.1831E-05*(x(2)*180/pi)^2 -...
        3.1917E-04*(x(2)*180/pi) + 1.6834E-04 +...
        m*g/(qb*S)*sin(x(3)) - (x(4))/(qb*S))^2 + ...
        (x(1) - x(3))^2);
    [x,fval,exitflag,output] = fminsearch(fnc, ...
                [0, 0, 0, 8000], optimset('TolX',1e-15));
    x_trim = [x(1) x(2) x(3) x(4)];
    disp('-------------------------------------')
    disp('----- ---- Trim function ----- -----')
    disp(sprintf('AoA           [°]: %5.3d',x_trim(1)*180/pi))
    disp(sprintf('Elevator      [°]: %5.3d',x_trim(2)*180/pi))
    disp(sprintf('Pitch angle   [°]: %5.3d',x_trim(3)*180/pi))
    disp(sprintf('Engine thrust [N]: %6.1f',x_trim(4)))
    disp('-------------------------------------')
end

% ulozeni velicin do mat souboru
save(cfgmatfile);
