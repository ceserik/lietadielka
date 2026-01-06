% FLYmodel.m 
% Initialization script 
% Pavel Hospodar, CVUT, FEL 2010

%clc
%close all
%clear all



% inputs parameters
plane = input_param;    

% Aerodynamic and mass model initialization 
% Aircraft triming 
x_trim = trimConfig(plane);

% initial states 
H_init = plane.h;
v0 = plane.v0;
V_init = [plane.v0*cos(x_trim(1)) 0 plane.v0*sin(x_trim(1))]';
alpha_init  = 1*x_trim(1);
el_init     = x_trim(2);
pitch_init  = 1*x_trim(3);
eng_init    = x_trim(4);
Ts = 1/50;
Tsim = plane.Tsim;

% longitudinal state-space model
[Ap,Bp,Cp,Dp] = longMotion(plane.v0,plane.h,plane.m,...
    plane.Jzz,plane.x_T,alpha_init,pitch_init);

% lateral state-space model
[As,Bs,Cs,Ds] = latMotion(plane.v0,plane.h,plane.m,plane.Jxx,...
    plane.Jyy,plane.Jxy,alpha_init,pitch_init);

%figure(1)
% scatter(real(eig(Ap)),imag(eig(Ap)))
% hold on

if ~exist("legends")
    legends(1,1) =plane.v0;
    legends(1,2) =plane.h;
else
    legends(end+1,1) =plane.v0;
    legends(end,2) =plane.h;
end
analysis







 