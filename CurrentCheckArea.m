function [curr_u,curr_v] = CurrentCheckArea(box,startlat,startlong,i,j,m,datapoint_u,datapoint_v,long_map,lat_map,curr_u,curr_v)

    if box == 3
        curr_u(3) = datapoint_u(i-1,j,m);
        curr_v(3) = datapoint_v(i-1,j,m);
    elseif box == 4
        curr_u(4) = datapoint_u(i-1,j+1,m);
        curr_v(4) = datapoint_v(i-1,j+1,m);
    elseif box == 5
        curr_u(5) = datapoint_u(i,j+1,m);
        curr_v(5) = datapoint_v(i,j+1,m);
    elseif box == 6
        curr_u(6) = datapoint_u(i+1,j+1,m);
        curr_v(6) = datapoint_v(i+1,j+1,m);
    elseif box == 7 
        curr_u(7) = datapoint_u(i+1,j,m);
        curr_v(7) = datapoint_v(i+1,j,m);
    elseif box == 8
        curr_u(8) = datapoint_u(i+1,j-1,m);
        curr_v(8) = datapoint_v(i+1,j-1,m);
    elseif box == 1
       curr_u(1) = datapoint_u(i,j-1,m);
       curr_v(1) = datapoint_v(i,j-1,m);
    elseif box == 2
       curr_u(2) = datapoint_u(i-1,j-1,m);
       curr_v(2) = datapoint_v(i-1,j-1,m); 
    end 
end



