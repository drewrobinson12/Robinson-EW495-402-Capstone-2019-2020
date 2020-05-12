function [wsidot_u,wsidot_v,wind_from,flag] = WindCheckArea(box,startlat,startlong,i,j,u_w_map,v_w_map,long_map,lat_map,wsidot_u,wsidot_v)
    
    flag = 0;
    if startlat == 64 && startlong < 10
        flag = 2;
        wsidot_u = [0,0,0,0,0,0,0,0,0];
        wsidot_v = [0,0,0,0,0,0,0,0,0];
        wind_from = 90;
        return
    elseif startlat == 64 && startlong >= 10
        flag = 1;
        wsidot_u = [0,0,0,0,0,0,0,0,0];
        wsidot_v = [0,0,0,0,0,0,0,0,0];
        wind_from = 90;
        return
    elseif startlat == 1
        flag = 1;
        wsidot_u = [0,0,0,0,0,0,0,0,0];
        wsidot_v = [0,0,0,0,0,0,0,0,0];
        wind_from = 90;
        return
    end
        
    if box == 3
        boat_true = HeadingAdjuster(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(startlat-1,startlong),long_map(startlat-1,startlong));
        wind_from = WindFrom(v_w_map(i-1,j),u_w_map(i-1,j));
        [app_wind,sf] = AppWind(boat_true,wind_from);
        mag_wind_spd = hypot(v_w_map(i-1,j),u_w_map(i-1,j));
        wsidot_u(3) = sf*mag_wind_spd*cos(boat_true);
        wsidot_v(3) = sf*mag_wind_spd*sin(boat_true);
    elseif box == 4
        boat_true = HeadingAdjuster(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(startlat-1,startlong+1),long_map(startlat-1,startlong+1));
        wind_from = WindFrom(v_w_map(i-1,j+1),u_w_map(i-1,j+1));
        [app_wind,sf] = AppWind(boat_true,wind_from);
        mag_wind_spd = hypot(v_w_map(i-1,j+1),u_w_map(i-1,j+1));
        wsidot_u(4) = sf*mag_wind_spd*cos(boat_true);
        wsidot_v(4) = sf*mag_wind_spd*sin(boat_true);
    elseif box == 5
        boat_true = HeadingAdjuster(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(startlat,startlong+1),long_map(startlat,startlong+1));
        wind_from = WindFrom(v_w_map(i,j+1),u_w_map(i,j+1));
        [app_wind,sf] = AppWind(boat_true,wind_from);
        mag_wind_spd = hypot(v_w_map(i,j+1),u_w_map(i,j+1));
        wsidot_u(5) = sf*mag_wind_spd*cos(boat_true);
        wsidot_v(5) = sf*mag_wind_spd*sin(boat_true);
    elseif box == 6
        boat_true = HeadingAdjuster(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(startlat+1,startlong+1),long_map(startlat+1,startlong+1));
        wind_from = WindFrom(v_w_map(i+1,j+1),u_w_map(i+1,j+1));
        [app_wind,sf] = AppWind(boat_true,wind_from);
        mag_wind_spd = hypot(v_w_map(i+1,j+1),u_w_map(i+1,j+1));
        wsidot_u(6) = sf*mag_wind_spd*cos(boat_true);
        wsidot_v(6) = sf*mag_wind_spd*sin(boat_true);
    elseif box == 7 
        boat_true = HeadingAdjuster(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(startlat+1,startlong),long_map(startlat+1,startlong));
        wind_from = WindFrom(v_w_map(i+1,j),u_w_map(i+1,j));
        [app_wind,sf] = AppWind(boat_true,wind_from);
         mag_wind_spd = hypot(v_w_map(i+1,j),u_w_map(i+1,j));
        wsidot_u(7) = sf*mag_wind_spd*cos(boat_true);
        wsidot_v(7) = sf*mag_wind_spd*sin(boat_true);
    elseif box == 8
        boat_true = HeadingAdjuster(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(startlat+1,startlong-1),long_map(startlat+1,startlong-1));
        wind_from = WindFrom(v_w_map(i+1,j-1),u_w_map(i+1,j-1));
        [app_wind,sf] = AppWind(boat_true,wind_from);
        mag_wind_spd = hypot(v_w_map(i+1,j-1),u_w_map(i+1,j-1));
        wsidot_u(8) = sf*mag_wind_spd*cos(boat_true);
        wsidot_v(8) = sf*mag_wind_spd*sin(boat_true);
    elseif box == 1
        boat_true = HeadingAdjuster(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(startlat,startlong-1),long_map(startlat,startlong-1));
        wind_from = WindFrom(v_w_map(i,j-1),u_w_map(i,j-1));
        [app_wind,sf] = AppWind(boat_true,wind_from);
        mag_wind_spd = hypot(v_w_map(i,j-1),u_w_map(i,j-1));
        wsidot_u(1) = sf*mag_wind_spd*cos(boat_true);
        wsidot_v(1) = sf*mag_wind_spd*sin(boat_true);
    elseif box == 2
        boat_true = HeadingAdjuster(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(startlat-1,startlong-1),long_map(startlat-1,startlong-1));
        wind_from = WindFrom(v_w_map(i-1,j-1),u_w_map(i-1,j-1));
        [app_wind,sf] = AppWind(boat_true,wind_from);
        mag_wind_spd = hypot(v_w_map(i-1,j-1),u_w_map(i-1,j-1));
        wsidot_u(2) = sf*mag_wind_spd*cos(boat_true);
        wsidot_v(2) = sf*mag_wind_spd*sin(boat_true);
    end 
end

