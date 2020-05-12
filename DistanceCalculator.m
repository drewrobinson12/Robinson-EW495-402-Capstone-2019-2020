function [dist] = DistanceCalculator(lat_now,long_now,lat_next,long_next)
    y_nm = (lat_next-lat_now)*60; %dist in NM from lat to next lat
    x_low = 0.0034987*(lat_now)^2 - 1.3302*(lat_now) + 96.413; %long dist in NM at given lat
    x_hi = 0.0034987*(lat_next)^2 - 1.3302*(lat_next) + 96.413; %dist in nm at next lat
    x_nm = (long_next-long_now)*((x_low + x_hi)/2); %dist in NM from given long to final long (average bw both points)
    dist = hypot(x_nm,y_nm); %straight line distance (NM)
end

