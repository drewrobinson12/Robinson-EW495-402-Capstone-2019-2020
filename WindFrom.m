function [true_wind_from] = WindFrom(v_wind,u_wind)
    wind_to = atan2d(v_wind,u_wind);
    
    if wind_to < 0
        wind_to = wind_to + 360;
    end 

    if wind_to >= 0 & wind_to <= 180
        true_wind_from = wind_to + 180;
    else
        true_wind_from = wind_to - 180;
    end 
end

