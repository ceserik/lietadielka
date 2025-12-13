function plane = input_param
% function for selecting the working point of the flight envelope
disp('-----------------------------------')
disp(sprintf(' Select poin from flight envelope?'))
%nn = input(' [0 - no] [1 - yes]:');
nn = 1
if nn == 1 
    err_init = 1;
    Obalka = [0 0.168; 1000 0.182; 2000 0.196; 3000 0.211; 4000 0.228 
        5000 0.246; 6000 0.266; 7000 0.288; 8000 0.313; 9000 0.342 
       10000 0.380; 11000 0.435; 12000 0.505; 12732 0.593; 12884 0.638 
       12884 0.638; 12732 0.730; 12000 0.762; 11000 0.777; 10000 0.784
        9000 0.788;  8000 0.790; 7000 0.790; 6000 0.788; 5000 0.787
        4000 0.784; 3000 0.781; 2000  0.777; 1000  0.770; 0  0.764];

    MaxSpeed = [0 0.76; 1150 0.82; 13000 0.82];

    figure(1)
    hold on
    plot(Obalka(:,2),Obalka(:,1),'color',[53 151 0]/256,'linewidth',2)
    plot(MaxSpeed(:,2),MaxSpeed(:,1),'k','linewidth',2)
    xlabel('Mach number [-]')
    ylabel('altitude [m]')
    title('Flight envelope')
    grid on
    axis square

    waitforbuttonpress
    u = get(gca,'currentpoint');
    u = u(1,[1 2])';

    h = ishold;
    plot(u(1),u(2),'rx','linewidth',2)
    hold off
    while err_init
        prompt={'Air speed (68 - 280) [m.s^{-1}]:',...
                'Altitude (0 - 13000) [m]:',...
                'Aircraft mass (4350 - 8000) [kg]:',...
                'Inertia moment Jxx (1000 - 10000) [kg.m^2]:',...
                'Inertia moment Jyy (30000 - 60000) [kg.m^2]:',...
                'Inertia moment Jxy (0 - 3000)   [kg.m^2]:',...
                'Inertia moment Jzz (31000 - 36000) [kg.m^2]:',...
                'Center of gravity (0.16 - 0.27) [kg.m^2]:',...
                'Simulation time [s]:'};
        name = 'Input parameters';
        numlines=1;
        
        % vypocet rychlsti zvuku (Machova cisla) na zaklade MSA
        vyska = u(2);
        T_0         = 15;          % teplota [°]
        dTdh        = 6.51122e-3;  % teplotni gradient [°/m]
        if vyska < 11000
           Tx     = T_0 - dTdh*vyska;
        else
           Tx     = -56.5;
        end
        ax = 331.3+(0.607*Tx);
        defaultanswer={num2str(round(u(1)*ax)),num2str(round(u(2))),'6500','5000',...
            '45000','2200','34000','0.245','50'};
        options.Resize='on';
        options.Interpreter='tex';
        answer=inputdlg(prompt,name,numlines,defaultanswer,options);

        % meze pro jednotlive veliciny 
        v0_min      = 68;       v0_max      = 280;
        height_min  = 0;        height_max  = +13000;
        m_min       = 4350;     m_max       = 8000;
        Jzz_min     = 31000;    Jzz_max     = 36000;
        Jxx_min     = 1e3;      Jxx_max     = 1e4;
        Jyy_min     = 3e4;      Jyy_max     = 6e4;
        Jxy_min     = 0;        Jxy_max     = 3e3;
        x_T_min     = 0.16;     x_T_max     = 0.27;

        v0 = str2double(answer{1});
        h = str2double(answer{2});
        m = str2double(answer{3});
        Jxx = str2double(answer{4});
        Jyy = str2double(answer{5});
        Jxy = str2double(answer{6});
        Jzz = str2double(answer{7});
        x_T = str2double(answer{8});
        Tsim =  str2double(answer{9});
        
        % Kontrola vstupnich hodnot 
        if ((v0_min > v0) || ( v0 > v0_max ))
            disp(sprintf('%s od: %g do: %g \n',...
                                   'Hodnota v0 mimo povolenou mez',v0_min,v0_max));
            err_init = err_init + 1;
        end
        if ((height_min > h) || ( h > height_max ))
            disp(sprintf('%s od: %g do: %g \n',...
                       'Hodnota height mimo povolenou mez',height_min,height_max));
            err_init = err_init + 1;
        end
        if ((m_min > m) || ( m > m_max ))
            disp(sprintf('%s od: %g do: %g \n',...
                                      'Hodnota m mimo povolenou mez',m_min,m_max));
            err_init = err_init + 1;
        end
        if ((Jzz_min > Jzz) || ( Jzz > Jzz_max ))
            disp(sprintf('%s od: %g do: %g \n',...
                                'Hodnota Jzz mimo povolenou mez',Jzz_min,Jzz_max));
            err_init = err_init + 1;
        end
        if ((x_T_min > x_T) || ( x_T > x_T_max ))
            disp(sprintf('%s od: %g do: %g \n', ...
                                'Hodnota x_T mimo povolenou mez',x_T_min,x_T_max));
            err_init = err_init + 1;
        end
        if ((Jxx_min > Jxx) || ( Jxx > Jxx_max ))
            disp(sprintf('%s od: %g do: %g \n',...
                                'Hodnota Jxx mimo povolenou mez',Jxx_min,Jxx_max));
            err_init = err_init + 1;
        end
        if ((Jyy_min > Jyy) || ( Jyy > Jyy_max ))
            disp(sprintf('%s od: %g do: %g \n',...
                                'Hodnota Jyy mimo povolenou mez',Jyy_min,Jyy_max));
            err_init = err_init + 1;
        end
        if ((Jxy_min > Jxy) || ( Jxy > Jxy_max ))
            disp(sprintf('%s od: %g do: %g \n',...
                                'Hodnota Jxy mimo povolenou mez',Jxy_min,Jxy_max));
            err_init = err_init + 1;
        end

        if err_init == 1
            err_init = 0;
        else 
            err_init = 1;
        end
    end
else
    answerX={'240','3950','6500','5000','45000','2200','34000','0.245','50'};
    
    v0 = str2double(answerX{1});
    h = str2double(answerX{2});
    m = str2double(answerX{3});
    Jxx = str2double(answerX{4});
    Jyy = str2double(answerX{5});
    Jxy = str2double(answerX{6});
    Jzz = str2double(answerX{7});
    x_T = str2double(answerX{8});
    Tsim =  str2double(answerX{9});        
end


vyska = h;
T_0         = 15;          % teplota [°]
dTdh        = 6.51122e-3;  % teplotni gradient [°/m]
rho_0       = 1.2250;      % hustota vzduchu [kg/m^3]
if vyska < 11000
	T     = T_0 - dTdh*vyska;
else
	T     = -56.5;
end
a = 331.3+(0.607*T);
rho = rho_0*(1-vyska/44308)^4.2553;

    
plane.m   = m;      % hmotnost letadla
plane.h   = h;      % vyska letu
plane.v0  = v0;     % rychlost letu
plane.Jxx = Jxx;    % moment setrvacnosti pro osu x
plane.Jyy = Jyy;    % moment setrvacnosti pro osu y
plane.Jzz = Jzz;    % moment setrvacnosti pro osu z
plane.Jxy = Jxy;    % moment setrvacnosti pro osu z
plane.x_T = x_T;    % centrum teziste
plane.M = v0/a;     % Machovo cislo
plane.Tsim = Tsim;  % delka simulace
plane.rho = rho;    % hustota vzduchu

disp('-----------------------------------')
disp(sprintf(' Air speed:     %5.1f',v0))
disp(sprintf(' Mach number:   %5.3f',v0/a))
disp(sprintf(' Altitude:      %5.1f',h))
disp(sprintf(' Aircraft mass: %5.1d',m))
disp(sprintf(' Inertia moment Jxx:%5.1d',Jxx))
disp(sprintf(' Inertia moment Jyy:%5.1d',Jyy))
disp(sprintf(' Inertia moment Jzz:%5.1d',Jzz))
disp(sprintf(' Inertia moment Jxy:%5.1d',Jxy))
disp(sprintf(' Center of gravity:%5.3g',x_T))
disp('-----------------------------------')
