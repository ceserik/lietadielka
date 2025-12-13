function [A,B,C,D] = longMotion(v0,h,m,Jzz,x_T,alpha0,theta0)
disp('-----------------------------------')
disp(' Dynamic model of longitudinal motion');


% Zakladni konstanty
g0          = 9.81; 
R_E         = 6378160;
T_0         = 288.15;
dTdh        = -6.51122e-3;
rho_0       = 1.2250;
S           = 18.8; 
b           = 2.15;
x_T0        = 0.245;
v0_min      = 68;
v0_max      = 280;
height_min  = -300;
height_max  = +13000;
m_min       = 4350;
m_max       = 8000;
Jzz_min     = 31000;
Jzz_max     = 36000;
x_T_min     = 0.16;
x_T_max     = 0.27;
x_cy0       = [0,0.4,0.6,0.7,0.75,0.8,0.85,0.9]; 
y_cy0       = [0.09,0.09,0.09,0.105,0.13,0.12,0.07,0.02];
y_mz0       = [0.025,0.025,0.03,0.025,0.02,0.05,-0.02,-0.03];
tx_M        = [0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.9];
tc_x_v      = [-0.174, -0.174, -0.0474, -0.0349, -0.0168, -0.0152, ...
               -0.075, -0.075];
tc_x_a      = [-4.44, -4.44, -2.44, -1.593, -1.215, 0.121, 0.106, ...
                0.106];
tc_x_d      = [0.025, 0.025, 0.018, 0.017, 0.018, 0.023, 0.023, ...
               0.023];
tc_y_v      = [0.652, 0.652, 0.0707, 0.079, 0.095, 0.1163, 0.1356, ...
               0.1356];
tc_y_ad     = [2.5, 2.5, 1.75, 1.25, 1.2, 0.65, 0.0, 0.0];
tc_y_w      = [4.45, 4.45, 4.5, 4.6, 4.7, 4.85, 4.75, 4.75];
tm_z_v      = [0.0042, 0.0042, 0.0042, 0.0057, 0.0, 0.0, 0.0, 0.0];
%
%--------------------------------------------------------------------------
% Kontrola vstupnich hodnot 
if ((v0_min > v0) || ( v0 > v0_max ))
    sprintf('%s od: %g do: %g \n',...
                            'Hodnota v0 mimo povolenou mez',v0_min,v0_max);
end
% 
if ((height_min > h) || ( h > height_max ))
    sprintf('%s od: %g do: %g \n',...
                'Hodnota height mimo povolenou mez',height_min,height_max);
end
% 
if ((m_min > m) || ( m > m_max ))
    sprintf('%s od: %g do: %g \n',...
                               'Hodnota m mimo povolenou mez',m_min,m_max);
end
%
if ((Jzz_min > Jzz) || ( Jzz > Jzz_max ))
    sprintf('%s od: %g do: %g \n',...
                         'Hodnota Jzz mimo povolenou mez',Jzz_min,Jzz_max);
end
%
if ((x_T_min > x_T) || ( x_T > x_T_max ))
    sprintf('%s od: %g do: %g \n', ...
                         'Hodnota x_T mimo povolenou mez',x_T_min,x_T_max);
end
% 
% -------------------------------------------------------------------------
% Vypoèty dle MSA
    T_0         = 15;          % teplota [°]
    dTdh        = 6.51122e-3;  % teplotni gradient [°/m]
    rho_0       = 1.2250;      % hustota vzduchu [kg/m^3]

    if h < 11000
       T     = T_0 - dTdh*h;
    else
       T     = -56.5;
    end

    rho_h = rho_0*(1-h/44308)^4.2553;
    a_h = 331.3+(0.607*T);
%     g = 9.81;
% -------------------------------------------------------------------------
% Vypocet alpha0 & theta0 
Mach_number = v0/a_h;
q           = 0.5*rho_h*(v0^2);
%
cy0 = interp1(x_cy0,y_cy0,Mach_number);
mz0 = interp1(x_cy0,y_mz0,Mach_number);
%
if (Mach_number <= 0.3)
    c_y_alpha = 4.4418;
else
    c_y_alpha = -62.716*(Mach_number^5)...
                 +122.09*(Mach_number^4) - 85.797*(Mach_number^3)...
                 +28.153*(Mach_number^2) - 3.8433*(Mach_number)...
                 + 4.541;
end
%
if (Mach_number <= 0.3)
    m_z_alpha = -0.4899;
else
    m_z_alpha = 5.9964*(Mach_number^4)...
                -10.271*(Mach_number^3) + 6.7921*(Mach_number^2)...
                -2.0056*(Mach_number) - 0.2708;
end

c_y_dv              = +0.305;
m_z_dv              = -0.88;
%
cy          = m*g0/q/S;
detA        = c_y_alpha*m_z_dv - c_y_dv*m_z_alpha;
detA_alpha  = (cy - cy0)*m_z_dv - (-mz0)*c_y_dv;
detA_dv     = c_y_alpha*(-mz0) - (cy - cy0)*m_z_alpha;

%--------------------------------------------------------------------------
% 6.Vypocet podelnych bezrozmernych aerodynamickych derivaci  
%
c_x_0			= 0.08;
c_y_0			= cy;
m_z_0			= mz0*0;
if (Mach_number <= 0.3)
	m_z_omega_z_bar	= -4.9934668;
	m_z_alpha_dot_bar   = -1.38271969;
else
	m_z_omega_z_bar	= +2.5926*(Mach_number^3)    - 8.7063*(Mach_number^2)... 
						  +4.8590*(Mach_number)    - 5.7376;
	m_z_alpha_dot_bar	= ...
        +502.427*(Mach_number^5)  - 1303.11*(Mach_number^4) ...
  	  +1330.2481*(Mach_number^3) - 668.0225*(Mach_number^2) ...
						  +164.08*(Mach_number)    - 17.0671;
end

c_x_v_bar			= interp1(tx_M,tc_x_v,Mach_number);
c_x_alpha			= -(10.081 - 23.378*Mach_number + 13.048*Mach_number^2);
c_x_omega_z_bar		= 0.0;
c_x_dv				= interp1(tx_M,tc_x_d,Mach_number);
c_y_v_bar			= interp1(tx_M,tc_y_v,Mach_number);
c_y_omega_z_bar		= interp1(tx_M,tc_y_w,Mach_number);
c_y_alpha_dot_bar	= interp1(tx_M,tc_y_ad,Mach_number);
m_z_v_bar   		= interp1(tx_M,tm_z_v,Mach_number); 

c_x_v_bar = 0;
c_y_v_bar = 0;
m_z_v_bar = 0;

% Coef.value = [cy0 c_y_alpha c_y_dv c_y_omega_z_bar c_y_alpha_dot_bar...
%     mz0 m_z_alpha m_z_dv m_z_omega_z_bar m_z_alpha_dot_bar ...
%     c_x_0 c_x_alpha c_x_dv]';
% Coef.name = {'c_y_0' 'c_y_alpha' 'c_y_dv' 'c_y_omega_z_bar' 'c_y_alpha_dot_bar'...
%     'm_z_0' 'm_z_alpha' 'm_z_dv' 'm_z_omega_z_bar' 'm_z_alpha_dot_bar' ...
%     'c_x_0' 'c_x_alpha' 'c_x_dv'}';
%-------------------------------------------------------------------------
% Vypocet rozmerovych aerodynamickych derivaci
tau				= 2*m/rho_h/S/v0;
X_v				= (c_x_v_bar + 2*c_x_0)*m/tau;
X_alpha			= (c_x_alpha + 2*c_x_0*alpha0)*v0*m/tau;
X_omega_z		= c_x_omega_z_bar*b*m/tau;
X_dv			= c_x_dv*v0*m/tau;
Y_v				= (c_y_v_bar + 2*c_y_0)*m/tau;
Y_alpha			= (c_y_alpha + 2*c_y_0*alpha0)*v0*m/tau;
Y_omega_z		= c_y_omega_z_bar*b*m/tau;
Y_alpha_dot		= c_y_alpha_dot_bar*b*m/tau;
Y_dv			= c_y_dv*v0*m/tau;
M_z_v			= (m_z_v_bar + 2*m_z_0)*m*b/tau;
M_z_alpha		= (m_z_alpha + 2*m_z_0*alpha0)*b*v0*m/tau;
M_z_omega_z     = m_z_omega_z_bar*(b^2)*m/tau;
M_z_alpha_dot   = m_z_alpha_dot_bar*(b^2)*m/tau;
M_z_dv			= m_z_dv*v0*b*m/tau;


%
%-------------------------------------------------------------------------
% Vypocet matic systemu
denom = m*v0 + Y_alpha_dot;
A = ...
[-X_v/m , -X_alpha/m , -(X_omega_z/m + alpha0*v0) , -g0*cos(theta0) ; ...
-Y_v/denom, -Y_alpha/denom, -(Y_omega_z - m*v0)/denom , ...
-m*g0*sin(theta0)/denom ; ...
(M_z_v + (x_T - x_T0)*b*Y_v - (M_z_alpha_dot + ...
(x_T - x_T0)*b*Y_alpha_dot)*Y_v/denom)/Jzz, ...
(M_z_alpha + (x_T - x_T0)*b*Y_alpha - (M_z_alpha_dot + ...
(x_T - x_T0)*b*Y_alpha_dot)*Y_alpha/denom)/Jzz, ...
(M_z_omega_z + (x_T - x_T0)*b*Y_omega_z - (M_z_alpha_dot ...
+ (x_T - x_T0)*b*Y_alpha_dot)*(Y_omega_z - m*v0)/denom)/Jzz, ...
-(M_z_alpha_dot + ...
(x_T - x_T0)*b*Y_alpha_dot)*m*g0*sin(theta0)/denom/Jzz; ...
0 , 0 , 1 , 0 ];

% % CHYBA V MATICI B
%     (M_z_dv + (x_T - x_T0)*b*Y_dv - (M_z_alpha_dot + ...
%     (x_T - x_T0)*b*Y_alpha_dot)/denom)/Jzz; 0];
% opraveno Pavel Hospodar
B = [-X_dv / m ; -Y_dv / denom ; ...
    (M_z_dv + (x_T - x_T0)*b*Y_dv - Y_dv*(M_z_alpha_dot + ...
    (x_T - x_T0)*b*Y_alpha_dot)/denom)/Jzz; 0];
C = [1 , 0 , 0 , 0; 
     0 , 1 , 0 , 0; ...
     0 , 0 , 1 , 0; ...
     0 , 0 , 0 , 1; ...
     -X_v / (m*g0) , -X_alpha / (m*g0) , -X_omega_z / (m*g0) , 0; ...
      Y_v / (m*g0) ,  Y_alpha / (m*g0) ,  Y_omega_z / (m*g0) , 0];
%
D = [ 0 ; 0 ; 0 ; 0 ; -X_dv / (m*g0) ; Y_dv / (m*g0) ];
%
%
%%
