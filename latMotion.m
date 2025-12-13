function [A,B,C,D] = latMotion(v0,height,m,Jxx,Jyy,Jxy,alpha0,theta0)
disp('-----------------------------------')
disp(' Dynamic model of lateral motion')
disp('-----------------------------------')
% Zakladni konstanty
R_E         = 6378160;
T_0         = 288.15;
dTdh        = -6.51122e-3;
rho_0       = 1.2250;
S           = 18.8; 
l           = 9.17;
S_rudder    = 0.7157;
b_rudder    = 0.4210;
L_SOP       = 4.0813;
k_SOP       = 0.95;
%
% -------------------------------------------------------------------------
% Vypocty dle MSA
    T_0         = 15;          % teplota [°]
    dTdh        = 6.51122e-3;  % teplotni gradient [°/m]
    rho_0       = 1.2250;      % hustota vzduchu [kg/m^3]
    vyska = height;
    if vyska < 11000
       T     = T_0 - dTdh*vyska;
    else
       T     = -56.5;
    end

    rho_h = rho_0*(1-vyska/44308)^4.2553;
    a_h = 331.3+(0.607*T);
    g0 = 9.81; 
    
    Mach_number = v0/a_h;
    q           = 1/2*rho_h*v0^2;
    cy = m*g0/q/S;

% -------------------------------------------------------------------------
% Vypocet bezrozmernych aerodynamickych derivaci
c_z_beta = -0.946 - 0.0333*cy + 0.0032*Mach_number - 0.16*Mach_number^2;
if Mach_number<0.45
    m_x_beta = -0.056 - 0.0576*cy - 0.0088*Mach_number;
else
    m_x_beta = -0.06 - 0.0576*cy - 0.0285*(Mach_number - 0.45);
end
m_y_beta = -0.173 + 0.01667*cy + 0.0092*Mach_number -0.1*Mach_number^2 +...
    0.0208*Mach_number^3;
m_x_omega_x = -0.5365 - 1.5267*Mach_number + 10.781*Mach_number^2 ...
    -20.666*Mach_number^3 + 12.4242*Mach_number^4;
if Mach_number<0.7
    m_y_omega_x = -0.07435 + 0.045*cy - 0.0092*Mach_number;
else
    m_y_omega_x = -0.08079 + 0.044*cy - 0.0333*(Mach_number-0.7);
end
m_x_omega_y = -0.075 - 0.215*cy -0.016*Mach_number;
m_y_omega_y = -0.2291 - 0.0065*cy -0.0236*cy^2 + 0.0069*cy^3 ...
    - 0.0219*Mach_number + 0.0097*Mach_number^2 - 0.0447*Mach_number^3;
c_z_dH = -0.1845 + 0.0325*cy -0.0093*Mach_number -0.0228*Mach_number^2 ...
    -0.012*Mach_number^3;
m_x_dH = -0.03 + 0.008*cy;
m_y_dH = 0.015*cy - 0.095 -0.0056*Mach_number -0.0089*Mach_number^2 ...
    - 0.0068*Mach_number^3;
c_z_dEE = -0.02;
if Mach_number<=0.6
    m_x_dEE = -0.156 + 0.04*(Mach_number -0.2);
else
    m_x_dEE = -1.1039 + 4.6696*Mach_number - 7.5941*Mach_number^2 ...
        + 4.1483*Mach_number^3;
end
m_y_dEE = -0.0015;

%-------------------------------------------------------------------------
% Vypocet rozmerovych aerodynamickych derivaci
Z_dEE       = 1*c_z_dEE*q*S;
Z_beta      = 1*c_z_beta*q*S;
Z_dH        = 1*c_z_dH*q*S;
M_x_beta    = 1*m_x_beta*q*S*l; 
M_y_beta    = 1*m_y_beta*q*S*l;
M_x_omega_x = 1*m_x_omega_x*l/2/v0*q*S*l; 
M_y_omega_x = 1*m_y_omega_x*l/2/v0*q*S*l;
M_x_omega_y = 1*m_x_omega_y*l/2/v0*q*S*l;
M_y_omega_y = 1*m_y_omega_y*l/2/v0*q*S*l;
M_x_dEE     = 1*m_x_dEE*q*S*l;
M_y_dEE     = 1*m_y_dEE*q*S*l;
M_x_dH      = 1*m_x_dH*q*S*l;
M_y_dH      = 1*m_y_dH*q*S*l; 

%-------------------------------------------------------------------------
% Transformace derivaci
denominator   = Jxx*Jyy - Jxy^2;
MM_x_beta     = Jyy/denominator*M_x_beta + Jxy/denominator*M_y_beta;
MM_y_beta     = Jxy/denominator*M_x_beta + Jxx/denominator*M_y_beta;
MM_x_omega_x  = Jyy/denominator*M_x_omega_x + Jxy/denominator*M_y_omega_x;
MM_y_omega_x  = Jxy/denominator*M_x_omega_x + Jxx/denominator*M_y_omega_x;
MM_x_omega_y  = Jyy/denominator*M_x_omega_y + Jxy/denominator*M_y_omega_y;
MM_y_omega_y  = Jxy/denominator*M_x_omega_y + Jxx/denominator*M_y_omega_y;
MM_x_dH       = Jyy/denominator*M_x_dH + Jxy/denominator*M_y_dH;
MM_y_dH       = Jxy/denominator*M_x_dH + Jxx/denominator*M_y_dH;
MM_x_dEE      = Jyy/denominator*M_x_dEE + Jxy/denominator*M_y_dEE;
MM_y_dEE      = Jxy/denominator*M_x_dEE + Jxx/denominator*M_y_dEE;

%-------------------------------------------------------------------------
% Vypocet matic systemu

A = ...
[Z_beta/m/v0 , tan(alpha0) , 1 , g0/v0*cos(theta0), 0;...
MM_x_beta,MM_x_omega_x, MM_x_omega_y, 0, 0;...
MM_y_beta,MM_y_omega_x, MM_y_omega_y, 0, 0;...
0        ,     1      , -tan(theta0), 0, 0;...
0        ,     0      ,  sec(theta0), 0, 0];
                                           
B = [Z_dEE/m/v0, Z_dH/m/v0; ...
     MM_x_dEE, MM_x_dH; ...
     MM_y_dEE, MM_y_dH; ...
     0, 0; ...
     0, 0];
 
C = [eye(5);Z_beta/m/g0 0 0 0 0];
D = [zeros(5,2);Z_dEE/m/g0 Z_dH/m/g0]; 


