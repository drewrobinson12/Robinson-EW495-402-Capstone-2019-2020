function [app_wind,sf] = AppWind(boat_true,wind_from_true)
        

        app_wind = wind_from_true - boat_true;
        if app_wind < 0 
            app_wind = app_wind + 360;
        end

        %{
        if app_wind > 30 & app_wind < 90
            %sf = 0.05;
        elseif (app_wind >= 330 & app_wind <= 360) || (app_wind >= 0 & app_wind <= 30)
           % sf = 0.0;
        elseif app_wind >= 270 & app_wind < 330
           % sf = 0.05;
        elseif app_wind >= 210 & app_wind < 270
           % sf = 0.075;
        elseif app_wind >= 150 & app_wind < 210
           % sf = 0.15;
        elseif app_wind >= 90 & app_wind < 150
           % sf = 0.075;
        end
        %}
        
        if (app_wind >= 335 && app_wind <= 360) || (app_wind >= 0 && app_wind < 25)
            sf = 0.0;
        elseif (app_wind >= 25 && app_wind < 40) || (app_wind >= 320 && app_wind < 335)
            sf = 0.1;
        elseif (app_wind >= 40 && app_wind < 60) || (app_wind >= 300 && app_wind < 320)
            sf = 0.116;
        elseif (app_wind >= 60 && app_wind < 80) || (app_wind >= 280 && app_wind < 300)
            sf = 0.123;
        elseif (app_wind >= 80 && app_wind < 100) || (app_wind >= 260 && app_wind < 280)
            sf = 0.121;
        elseif (app_wind >= 100 && app_wind < 120) || (app_wind >= 240 && app_wind < 260)
            sf = 0.116;
        elseif (app_wind >= 120 && app_wind < 140) || (app_wind >= 220 && app_wind < 240)
            sf = 0.116;
        elseif (app_wind >= 140 && app_wind < 160) || (app_wind >= 200 && app_wind < 220)
            sf = 0.119;  
        elseif (app_wind >= 160 && app_wind < 200)  
            sf = 0.123;
        end
        
end

