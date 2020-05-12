function [lat_path,long_path,seg_dist,total_dist,path_spd,seg_time,total_time,i,j,startlat,startlong,min_count,map_time,flag] = LeadIf2(k,m,i,j,lat_range,endlat_s,endlat_n,long_range,endlong_e,endlong_w,uv_map_lead,curr_angle,startlat,startlong,lat_map,long_map,total_dist,total_time,lat_path,long_path,b2t_lat,b2t_long,tgt_idx,datapoint_uw,datapoint_vw,datapoint_u,datapoint_v,min_count,seg_dist,seg_time,map_time)

result_u = [];
result_v = [];
result_mag = [];
result_dir = [];
wsidot_u = [];
wsidot_v = [];
curr_u = [];
curr_v = [];
spd_sf = 1.0;
flag = 0;
    if startlat == 64 && startlong < 10
        flag = 2;
        path_spd = 1;
        return
    elseif startlat == 64 && startlong >= 10
        flag = 1;
        path_spd = 1;
        return
    elseif startlat == 1
        flag = 1;
        path_spd = 1;
        return
    end

%*****BOX INFO COLLECTION*****
if j == 1 && i == find(lat_range==endlat_s) %condition for bottom left corner
    surr_curr_mag_lead = rmmissing([uv_map_lead(i-1,j),uv_map_lead(i-1,j+1),uv_map_lead(i,j+1)]); %top,upper right, right
    
    surr_curr_dir_lead = [curr_angle(i-1,j),curr_angle(i-1,j+1),curr_angle(i,j+1)];
    zind = find(surr_curr_dir_lead~=0);
    surr_curr_dir_lead = surr_curr_dir_lead(zind);

    for boxy = 3:5
        [wsidot_u,wsidot_v,wind_from,flag] = WindCheckArea(boxy,startlat,startlong,i,j,datapoint_uw,datapoint_vw,long_map,lat_map,wsidot_u,wsidot_v);
        [curr_u,curr_v] = CurrentCheckArea(boxy,startlat,startlong,i,j,m,datapoint_u,datapoint_v,long_map,lat_map,curr_u,curr_v);
        result_u(boxy) = curr_u(boxy) + wsidot_u(boxy);
        result_v(boxy) = curr_v(boxy) + wsidot_v(boxy);
        result_mag(boxy) = hypot(result_u(boxy),result_v(boxy));
        result_dir(boxy) = unit2comp(atan2(result_v(boxy),result_u(boxy)));
    end         
    mean_result_mag = mean(rmmissing(result_mag(3:5)));
    mean_result_dir = mean(rmmissing(result_dir(3:5)));
elseif j == 1 && i == find(lat_range==endlat_n) %condition for top left corner
    surr_curr_mag_lead = rmmissing([uv_map_lead(i,j+1),uv_map_lead(i+1,j+1),uv_map_lead(i+1,j)]); %right, bottom right, bottom
    
    surr_curr_dir_lead = [curr_angle(i,j+1),curr_angle(i+1,j+1),curr_angle(i+1,j)];
    zind = find(surr_curr_dir_lead~=0);
    surr_curr_dir_lead = surr_curr_dir_lead(zind);

    for boxy = 5:7
        [wsidot_u,wsidot_v,wind_from,flag] = WindCheckArea(boxy,startlat,startlong,i,j,datapoint_uw,datapoint_vw,long_map,lat_map,wsidot_u,wsidot_v);
        [curr_u,curr_v] = CurrentCheckArea(boxy,startlat,startlong,i,j,m,datapoint_u,datapoint_v,long_map,lat_map,curr_u,curr_v);
        result_u(boxy) = curr_u(boxy) + wsidot_u(boxy);
        result_v(boxy) = curr_v(boxy) + wsidot_v(boxy);
        result_mag(boxy) = hypot(result_u(boxy),result_v(boxy));
        result_dir(boxy) = unit2comp(atan2(result_v(boxy),result_u(boxy)));
    end        
    mean_result_mag = mean(rmmissing(result_mag(5:7)));
    mean_result_dir = mean(rmmissing(result_dir(5:7)));
elseif j == 1 && i~=find(lat_range==endlat_s) && i~=find(lat_range==endlat_n) %condition for left edge
    surr_curr_mag_lead = rmmissing([uv_map_lead(i-1,j),uv_map_lead(i-1,j+1),uv_map_lead(i,j+1),... %top, top right, right
        uv_map_lead(i+1,j+1),uv_map_lead(i+1,j)]); %bottom right, bottom
    
    surr_curr_dir_lead = [curr_angle(i-1,j),curr_angle(i-1,j+1),curr_angle(i,j+1),...
        curr_angle(i+1,j+1),curr_angle(i+1,j)];
    zind = find(surr_curr_dir_lead~=0);
    surr_curr_dir_lead = surr_curr_dir_lead(zind);

    for boxy = 3:7
        [wsidot_u,wsidot_v,wind_from,flag] = WindCheckArea(boxy,startlat,startlong,i,j,datapoint_uw,datapoint_vw,long_map,lat_map,wsidot_u,wsidot_v);
        [curr_u,curr_v] = CurrentCheckArea(boxy,startlat,startlong,i,j,m,datapoint_u,datapoint_v,long_map,lat_map,curr_u,curr_v);
        result_u(boxy) = curr_u(boxy) + wsidot_u(boxy);
        result_v(boxy) = curr_v(boxy) + wsidot_v(boxy);
        result_mag(boxy) = hypot(result_u(boxy),result_v(boxy));
        result_dir(boxy) = unit2comp(atan2(result_v(boxy),result_u(boxy)));
    end        
    mean_result_mag = mean(rmmissing(result_mag(3:7)));
    mean_result_dir = mean(rmmissing(result_dir(3:7)));
elseif i == 1 && j~=find(long_range==endlong_w) && j~= find(long_range==endlong_e) %condition for top edge
    surr_curr_mag_lead = rmmissing([uv_map_lead(i,j-1),uv_map_lead(i,j+1),... %left, right
        uv_map_lead(i+1,j+1),uv_map_lead(i+1,j),uv_map_lead(i+1,j-1)]); %bottom right, bottom, bottom left
    
    surr_curr_dir_lead = [curr_angle(i,j-1),curr_angle(i,j+1),... %left, right
        curr_angle(i+1,j+1),curr_angle(i+1,j),curr_angle(i+1,j-1)]; %bottom right, bottom, bottom left
    zind = find(surr_curr_dir_lead~=0);
    surr_curr_dir_lead = surr_curr_dir_lead(zind);

    for boxy = [5,6,7,8,1]
        [wsidot_u,wsidot_v,wind_from,flag] = WindCheckArea(boxy,startlat,startlong,i,j,datapoint_uw,datapoint_vw,long_map,lat_map,wsidot_u,wsidot_v);
        [curr_u,curr_v] = CurrentCheckArea(boxy,startlat,startlong,i,j,m,datapoint_u,datapoint_v,long_map,lat_map,curr_u,curr_v);
        result_u(boxy) = curr_u(boxy) + wsidot_u(boxy);
        result_v(boxy) = curr_v(boxy) + wsidot_v(boxy);
        result_mag(boxy) = hypot(result_u(boxy),result_v(boxy));
        result_dir(boxy) = unit2comp(atan2(result_v(boxy),result_u(boxy)));
    end        
    mean_result_mag = mean(rmmissing(result_mag([5,6,7,8,1])));
    mean_result_dir = mean(rmmissing(result_dir([5,6,7,8,1])));
elseif i == 1 && j == find(long_range==endlong_e) %condition for top right corner
    surr_curr_mag_lead = rmmissing([uv_map_lead(i,j-1),... %left
        uv_map_lead(i+1,j),uv_map_lead(i+1,j-1)]); %bottom, bottom left
    
    surr_curr_dir_lead = [curr_angle(i,j-1),... %left
        curr_angle(i+1,j),curr_angle(i+1,j-1)]; %bottom, bottom left
    zind = find(surr_curr_dir_lead~=0);
    surr_curr_dir_lead = surr_curr_dir_lead(zind);

    for boxy = [1,7,8]
        [wsidot_u,wsidot_v,wind_from,flag] = WindCheckArea(boxy,startlat,startlong,i,j,datapoint_uw,datapoint_vw,long_map,lat_map,wsidot_u,wsidot_v);
        [curr_u,curr_v] = CurrentCheckArea(boxy,startlat,startlong,i,j,m,datapoint_u,datapoint_v,long_map,lat_map,curr_u,curr_v);
        result_u(boxy) = curr_u(boxy) + wsidot_u(boxy);
        result_v(boxy) = curr_v(boxy) + wsidot_v(boxy);
        result_mag(boxy) = hypot(result_u(boxy),result_v(boxy));
        result_dir(boxy) = unit2comp(atan2(result_v(boxy),result_u(boxy)));
    end        
    mean_result_mag = mean(rmmissing(result_mag([1,7,8])));
    mean_result_dir = mean(rmmissing(result_dir([1,7,8])));
elseif i~=1 && i~=find(lat_range==endlat_s) && j == find(long_range==endlong_e) %condition for right edge
    surr_curr_mag_lead = rmmissing([uv_map_lead(i,j-1),uv_map_lead(i-1,j-1),uv_map_lead(i-1,j),... %left through top
        uv_map_lead(i+1,j),uv_map_lead(i+1,j-1)]); %bottom, bottom left
    
    surr_curr_dir_lead = [curr_angle(i,j-1),curr_angle(i-1,j-1),curr_angle(i-1,j),... %left through top
        curr_angle(i+1,j),curr_angle(i+1,j-1)]; %bottom, bottom left
    zind = find(surr_curr_dir_lead~=0);
    surr_curr_dir_lead = surr_curr_dir_lead(zind);

    for boxy = [6,7,8,1,2]
        [wsidot_u,wsidot_v,wind_from,flag] = WindCheckArea(boxy,startlat,startlong,i,j,datapoint_uw,datapoint_vw,long_map,lat_map,wsidot_u,wsidot_v);
        [curr_u,curr_v] = CurrentCheckArea(boxy,startlat,startlong,i,j,m,datapoint_u,datapoint_v,long_map,lat_map,curr_u,curr_v);
        result_u(boxy) = curr_u(boxy) + wsidot_u(boxy);
        result_v(boxy) = curr_v(boxy) + wsidot_v(boxy);
        result_mag(boxy) = hypot(result_u(boxy),result_v(boxy));
        result_dir(boxy) = unit2comp(atan2(result_v(boxy),result_u(boxy)));
    end        
    mean_result_mag = mean(rmmissing(result_mag([6,7,8,1,2])));
    mean_result_dir = mean(rmmissing(result_dir([6,7,8,1,2])));
elseif i == find(lat_range==endlat_s) && j == find(long_range==endlong_e) %condition for bottom right corner
    surr_curr_mag_lead = rmmissing([uv_map_lead(i,j-1),uv_map_lead(i-1,j-1),uv_map_lead(i-1,j)]); %left through top
    surr_curr_dir_lead = [curr_angle(i,j-1),curr_angle(i-1,j-1),curr_angle(i-1,j)]; %left through top
    zind = find(surr_curr_dir_lead~=0);
    surr_curr_dir_lead = surr_curr_dir_lead(zind);

    for boxy = 1:3
        [wsidot_u,wsidot_v,wind_from,flag] = WindCheckArea(boxy,startlat,startlong,i,j,datapoint_uw,datapoint_vw,long_map,lat_map,wsidot_u,wsidot_v);
        [curr_u,curr_v] = CurrentCheckArea(boxy,startlat,startlong,i,j,m,datapoint_u,datapoint_v,long_map,lat_map,curr_u,curr_v);
        result_u(boxy) = curr_u(boxy) + wsidot_u(boxy);
        result_v(boxy) = curr_v(boxy) + wsidot_v(boxy);
        result_mag(boxy) = hypot(result_u(boxy),result_v(boxy));
        result_dir(boxy) = unit2comp(atan2(result_v(boxy),result_u(boxy)));
    end        
    mean_result_mag = mean(rmmissing(result_mag(1:3)));
    mean_result_dir = mean(rmmissing(result_dir(1:3)));
elseif i == find(lat_range==endlat_s) && j~=find(long_range==endlong_w) && j~=find(long_range==endlong_e) %condition for bottom edge
    surr_curr_mag_lead = rmmissing([uv_map_lead(i,j-1),uv_map_lead(i-1,j-1),uv_map_lead(i-1,j),uv_map_lead(i-1,j+1),uv_map_lead(i,j+1)]); %left through right
    
    surr_curr_dir_lead = [curr_angle(i,j-1),curr_angle(i-1,j-1),curr_angle(i-1,j),curr_angle(i-1,j+1),curr_angle(i,j+1)]; %left through right
    zind = find(surr_curr_dir_lead~=0);
    surr_curr_dir_lead = surr_curr_dir_lead(zind);

    for boxy = 1:5
        [wsidot_u,wsidot_v,wind_from,flag] = WindCheckArea(boxy,startlat,startlong,i,j,datapoint_uw,datapoint_vw,long_map,lat_map,wsidot_u,wsidot_v);
        [curr_u,curr_v] = CurrentCheckArea(boxy,startlat,startlong,i,j,m,datapoint_u,datapoint_v,long_map,lat_map,curr_u,curr_v);
        result_u(boxy) = curr_u(boxy) + wsidot_u(boxy);
        result_v(boxy) = curr_v(boxy) + wsidot_v(boxy);
        result_mag(boxy) = hypot(result_u(boxy),result_v(boxy));
        result_dir(boxy) = unit2comp(atan2(result_v(boxy),result_u(boxy)));
    end    
    mean_result_mag = mean(rmmissing(result_mag(1:5)));
    mean_result_dir = mean(rmmissing(result_dir(1:5)));
else %regular, non-edge/corner decision
    surr_curr_mag_lead = rmmissing([uv_map_lead(i,j-1),uv_map_lead(i-1,j-1),uv_map_lead(i-1,j),uv_map_lead(i-1,j+1),uv_map_lead(i,j+1),... %left through right
        uv_map_lead(i+1,j+1),uv_map_lead(i+1,j),uv_map_lead(i+1,j-1)]); %bottom right, bottom, bottom left
    
    surr_curr_dir_lead = [curr_angle(i,j-1),curr_angle(i-1,j-1),curr_angle(i-1,j),curr_angle(i-1,j+1),curr_angle(i,j+1),... %left through right
        curr_angle(i+1,j+1),curr_angle(i+1,j),curr_angle(i+1,j-1)]; %bottom right, bottom, bottom left
    
    zind = find(surr_curr_dir_lead~=0);
    surr_curr_dir_lead = surr_curr_dir_lead(zind);

    for boxy = 1:8
        [wsidot_u,wsidot_v,wind_from,flag] = WindCheckArea(boxy,startlat,startlong,i,j,datapoint_uw,datapoint_vw,long_map,lat_map,wsidot_u,wsidot_v);
        [curr_u,curr_v] = CurrentCheckArea(boxy,startlat,startlong,i,j,m,datapoint_u,datapoint_v,long_map,lat_map,curr_u,curr_v);
        result_u(boxy) = curr_u(boxy) + wsidot_u(boxy);
        result_v(boxy) = curr_v(boxy) + wsidot_v(boxy);
        result_mag(boxy) = hypot(result_u(boxy),result_v(boxy));
        result_dir(boxy) = unit2comp(atan2(result_v(boxy),result_u(boxy)));
    end  
    mean_result_mag = mean(rmmissing(result_mag(1:8)));
    mean_result_dir = mean(rmmissing(result_dir(1:8)));
end

if flag == 2
    path_spd = 0;
    return
end

mean_surr_curr_mag_lead = mean(surr_curr_mag_lead); %kts
mean_surr_curr_dir_lead = mean(surr_curr_dir_lead); %degrees true

thresh = 1.1;
curr_sweep = 1.6;

if mean_surr_curr_mag_lead < thresh %THRESHOLD FOR MIN CURRENT
    min_count = min_count+1;
else
    min_count = 0;
end

%%(surr_curr_mag_lead)
%%(surr_curr_dir_lead)

% *****DECISION POINT*******

% if (mean_surr_curr_dir >= 337.5 && mean_surr_curr_dir <= 360) || (mean_surr_curr_dir >= 0 && mean_surr_curr_dir < 22.5)
box_bearings = [15,60,105,150];
check_stash = [];
check_stash2 = [];
if min_count <= 3 %IF 3 OR FEWER ITERATIONS IN A ROW HAVE CURRENT < MIN THRESHOLD
   %Remain in dominant current mode
    if mean_surr_curr_mag_lead > curr_sweep %SWEPT WITH CURRENT CONDITION
        for q = 3:6
            check_stash(q-2) = abs(box_bearings(q-2)-surr_curr_dir_lead(q)); %find how well current direction matches each boxy in 3-6
            [val,next_box] = min(check_stash); %find boxy with closest match to current direction
            check_stash2(q-2) = abs(box_bearings(q-2)-mean_surr_curr_dir_lead); %look at box compared to average
            [val2,next_box2] = min(check_stash2); %compare individual comparison and average comparison
                if val <= val2 %if one is closer fit than other
                    next_box = next_box; %choose closer fit of two
                else
                    next_box = next_box2; %choose closer fit of two
                end
            %('Swept w/ current > 1.0')
                 if mean_surr_curr_dir_lead > 315 && mean_surr_curr_dir_lead <= 360
                     next_box = 1;
                 end
                 path_spd = spd_sf*uv_map_lead(i,j);
            continue
        end
    elseif mean_surr_curr_mag_lead >= thresh %WIND ASSIST WITHIN DOMINANT CURRENT
        %IF WIND ASSIST MAKES TOTAL SPD GREATER AND IN DIRECTION OF DESIRED TRAVEL
        if mean_result_mag > curr_sweep && (mean_result_dir > 337.5 && mean_result_dir <= 360) || (mean_result_dir >= 0 && mean_result_dir <= 157.5)
            for q = 3:6
                %FIND DIFFERENCE BETWEEN BOX DIRECTION AND COMBINED DIRECTION
                check_stash(q-2) = abs(box_bearings(q-2)-result_dir(q)); 
                %GO TOWARDS BOX OF MINIMUM DIFFERENCE
                [val,next_box] = min(check_stash); %find boxy with lowest difference
                check_stash2(q-2) = abs(box_bearings(q-2)-mean_result_dir); %look at box compared to average
                [val2,next_box2] = min(check_stash2); %compare individual comparison and average comparison
                    if val <= val2 %if one is closer fit than other
                        next_box = next_box; %choose closer fit of two
                    else
                        next_box = next_box2; %choose closer fit of two
                    end
                    if mean_result_dir > 315 && mean_result_dir <= 360
                     next_box = 1;
                    end
                    path_spd = spd_sf*mean_result_mag;
                %('Wind assist in dominant current')
                continue
            end
        else %IF IT DOESN'T ASSIST WELL ENOUGH
            %RESORT TO SWEPT BY CURRENT CONDITION
            for q = 3:6
                check_stash(q-2) = abs(box_bearings(q-2)-surr_curr_dir_lead(q)); %find how well current direction matches each boxy in 3-6
                [val,next_box] = min(check_stash); %find boxy with closest match to current direction
                %('Swept w/ current bc bad assist')
                check_stash2(q-2) = abs(box_bearings(q-2)-mean_surr_curr_dir_lead); %look at box compared to average
                [val2,next_box2] = min(check_stash2); %compare individual comparison and average comparison
                    if val <= val2 %if one is closer fit than other
                        next_box = next_box; %choose closer fit of two
                    else
                        next_box = next_box2; %choose closer fit of two
                    end
                    path_spd = spd_sf*uv_map_lead(i,j);
                continue
            end
        end
        
    %WIND ASSIST WITHIN WEAK CURRENT
    elseif mean_surr_curr_mag_lead < thresh && mean_result_mag > thresh 
        for q = 3:6
            %FIND DIFFERENCE BETWEEN BOX DIRECTION AND COMBINED DIRECTION
            check_stash(q-2) = abs(box_bearings(q-2)-result_dir(q)); 
            %GO TOWARDS BOX OF MINIMUM DIFFERENCE
            [val,next_box] = min(check_stash); %find boxy with lowest difference
            %('Wind assist in weak current')
            check_stash2(q-2) = abs(box_bearings(q-2)-mean_result_dir); %look at box compared to average
            [val2,next_box2] = min(check_stash2); %compare individual comparison and average comparison
                if val <= val2 %if one is closer fit than other
                    next_box = next_box; %choose closer fit of two
                else
                    next_box = next_box2; %choose closer fit of two
                end
                path_spd = spd_sf*mean_result_mag;
            continue
        end
        %}
    
    %IF CURRENT AND WIND ASSIST ARE BOTH WEAK
    elseif (mean_surr_curr_dir_lead > 315 && mean_surr_curr_dir_lead <= 360) || (mean_surr_curr_dir_lead >= 0 && mean_surr_curr_dir_lead <= 157.5) %if current pushing to boxes 1,2,3, or 4
        %CALCULATE BEARING TO TARGET
        bearing2tgt = HeadingAdjuster(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(b2t_lat(tgt_idx),1),long_map(1,b2t_long(tgt_idx))); %calc bearing to finish line
        %MAY WANT TO CHECK STASH WITH COMBINED AND THEN CHOOSE BETTER OF 2
        for q = 3:6
            check_stash(q-2) = abs(bearing2tgt-surr_curr_dir_lead(q)); %find difference bw bearing to finish and current direction for each boxy
            [val,next_box] = min(check_stash); %find boxy with lowest difference
            check_stash2(q-2) = abs(bearing2tgt-mean_surr_curr_dir_lead); %look at box compared to average
            [val2,next_box2] = min(check_stash2); %compare individual comparison and average comparison
                if val <= val2 %if one is closer fit than other
                    next_box = next_box; %choose closer fit of two
                else
                    next_box = next_box2; %choose closer fit of two
                end
            %('Both weak, but 1-4')
            path_spd = spd_sf*mean_result_mag;
        end
    %IF CURRENT PUSHING SOUTH
    elseif mean_surr_curr_dir_lead > 157.5 && mean_surr_curr_dir_lead <= 202.5
        %SAIL EAST
        next_box = 3;
        path_spd = spd_sf*mean_result_mag;
        %('Weak pushing south')
    %IF CURRENT PUSHING SOUTHWEST
    elseif mean_surr_curr_dir_lead > 202.5 && mean_surr_curr_dir_lead <= 247.5 %if current pushing to boxy 8
        %SAIL SOUTHEAST
        next_box = 4; %sail perpendicular
        path_spd = spd_sf*mean_result_mag;
        %('Weak pushing southwest')
    %IF CURRENT PUSHING WEST
    elseif mean_surr_curr_dir_lead > 247.5 && mean_surr_curr_dir_lead <= 292.5 %if current pushing to boxy 1
        %SAIL NORTH
        next_box = 1; %sail perpendicular
        path_spd = spd_sf*mean_result_mag;
        %('Weak pushing north')
    %IF CURRENT PUSHING NORTHWEST
    elseif mean_surr_curr_dir_lead > 292.5 && mean_surr_curr_dir_lead <= 337.5 %if current pushing to boxy 2
        %SAIL NORTHEAST
        next_box = 2; %sail perpendicular
        path_spd = spd_sf*mean_result_mag;
        %('Weak pushing northwest')
    else
        %('Error')
        path_spd = 1;
        flag = 1;
        return
    end
    
%***IF MIN COUNT GREATER THAN THRESHOLD
%ENTER COMBINED MODE****
else
    check_stash = [];
    check_stash2 = [];
    if (mean_result_dir > 315 && mean_result_dir <= 360) || (mean_result_dir >= 0 && mean_result_dir <= 157.5) %if current pushing to boxes 1,2,3, or 4
    bearing2tgt = HeadingAdjuster(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(b2t_lat(tgt_idx),1),long_map(1,b2t_long(tgt_idx))); %calc bearing to finish line
    %MAY WANT TO CHECK STASH WITH COMBINED AND THEN CHOOSE BETTER OF 2
    for q = 3:6
        check_stash(q-2) = abs(bearing2tgt-result_dir(q)); %find difference bw bearing to finish and current direction for each boxy
        [val,next_box] = min(check_stash); %find boxy with lowest difference
        check_stash2(q-2) = abs(bearing2tgt-mean_result_dir); %look at box compared to average
        [val2,next_box2] = min(check_stash2); %compare individual comparison and average comparison
            if val <= val2 %if one is closer fit than other
                next_box = next_box; %choose closer fit of two
            else
                next_box = next_box2; %choose closer fit of two
            end
            path_spd = spd_sf*mean_result_mag;
        %('Counter on--both weak, 1-4')
        continue
    end
    elseif mean_result_dir > 157.5 && mean_result_dir <= 202.5
        %SAIL EAST
        next_box = 3;
        path_spd = spd_sf*mean_result_mag;
        %('Counter on--Weak pushing south')
    %IF CURRENT PUSHING SOUTHWEST
    elseif mean_result_dir > 180 && mean_result_dir <= 270 %if current pushing to boxy 8
        b2t_4 = HeadingAdjuster(lat_map(startlat+1,startlong+1),long_map(startlat+1,startlong+1),lat_map(b2t_lat(tgt_idx),1),long_map(1,b2t_long(tgt_idx))); %calc bearing to finish line
        b2t_1 = HeadingAdjuster(lat_map(startlat-1,startlong),long_map(startlat-1,startlong),lat_map(b2t_lat(tgt_idx),1),long_map(1,b2t_long(tgt_idx))); %calc bearing to finish line
        d2t_4 = DistanceCalculator(lat_map(startlat+1,startlong+1),long_map(startlat+1,startlong+1),lat_map(b2t_lat(tgt_idx),1),long_map(1,b2t_long(tgt_idx)));
        d2t_1 = DistanceCalculator(lat_map(startlat-1,startlong),long_map(startlat-1,startlong),lat_map(b2t_lat(tgt_idx),1),long_map(1,b2t_long(tgt_idx)));
        d2t_3 = DistanceCalculator(lat_map(startlat,startlong+1),long_map(startlat,startlong+1),lat_map(b2t_lat(tgt_idx),1),long_map(1,b2t_long(tgt_idx)));
        %IF BETTER BEARING TO TARGET IS TO TACK NORTHWEST THAN GO SOUTHEAST
        if d2t_1 < d2t_4 %&& d2t_1 < d2t_3
            next_box = 1; %1 close enough to 2 and avoids endless while
            %('Counter on--1 by decision')
        elseif d2t_4 < d2t_1 %&& d2t_4 < d2t_3
            next_box = 4; %sail perpendicular
             %('Counter on--4 by decision')
       % elseif d2t_3 <= d2t_1 && d2t_3 <= d2t_4
           % next_box = 3;
        end
        path_spd = spd_sf*mean_result_mag;
        
    %IF CURRENT PUSHING WEST
    elseif mean_result_dir > 247.5 && mean_result_dir <= 292.5 %if current pushing to boxy 1
        d2t_4 = DistanceCalculator(lat_map(startlat+1,startlong+1),long_map(startlat+1,startlong+1),lat_map(b2t_lat(tgt_idx),1),long_map(1,b2t_long(tgt_idx)));
        d2t_1 = DistanceCalculator(lat_map(startlat-1,startlong),long_map(startlat-1,startlong),lat_map(b2t_lat(tgt_idx),1),long_map(1,b2t_long(tgt_idx)));
        if d2t_1 < d2t_4 %&& d2t_1 < d2t_3
            next_box = 1; %1 close enough to 2 and avoids endless while
            %('Counter on--1 by decision')
        elseif d2t_4 < d2t_1 %&& d2t_4 < d2t_3
            next_box = 4; %sail perpendicular
             %('Counter on--4 by decision')
        %elseif d2t_3 <= d2t_1 && d2t_3 <= d2t_4
           % next_box = 3;
        end
        %SAIL NORTH
        %next_box = 1; %sail perpendicular
        path_spd = spd_sf*mean_result_mag;
        %('Counter on--Weak pushing north')
    %IF CURRENT PUSHING NORTHWEST
    elseif mean_result_dir > 292.5 && mean_result_dir <= 337.5 %if current pushing to boxy 2
        %SAIL NORTHEAST
        next_box = 2; %sail perpendicular
        path_spd = spd_sf*mean_result_mag;
        %('Counter on--Weak pushing northwest')
    else
        %('Error')
        path_spd = 1;
        flag = 1;
        return
        %next_box = 3; %default==sail east
    end
end

% %(check_stash)
% %(check_stash2)
%(next_box)
% %(next_box2)
%(strcat('Wind From: ',num2str(wind_from)))
%(strcat('Mag WSIDOT: ',num2str(hypot(wsidot_u,wsidot_v))))
%(strcat('Dir WSIDOT: ', num2str(unit2comp(atan2(wsidot_v,wsidot_u)))))
%(strcat('Mean combined mag: ',num2str(mean_result_mag)))
%(strcat('Mean combined dir: ',num2str(mean_result_dir)))


%****NEXT BOX "DO WHAT"*****

if next_box == 1
    max_index = 1;
    long_path(m,k) = long_map(startlat-1,startlong); %long of desired gradient
    lat_path(m,k) = lat_map(startlat-1,startlong); %lat of desired gradient
    seg_dist = SegmentDistance(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(startlat-1,startlong),long_map(startlat-1,startlong)); %NM
    total_dist = total_dist + seg_dist; %NM
    %path_spd = (uv_map(i,j)+uv_map(i-1,j))/2; %point to point speed, avg of start and endpoint (kts)
    %path_spd = uv_map_lead(i,j); %simplified to get around index error
    seg_time = seg_dist/path_spd; %hours
    map_time = map_time + seg_time; %hours
    total_time = total_time + seg_time; %hours
    i = i-1;
    startlat = startlat-1;
    j = j;
    startlong = startlong;
elseif next_box == 2
    max_index = 2;
    long_path(m,k) = long_map(startlat-1,startlong+1); %long of desired gradient
    lat_path(m,k) = lat_map(startlat-1,startlong+1); %lat of desired gradient
    seg_dist = SegmentDistance(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(startlat-1,startlong+1),long_map(startlat-1,startlong+1));
    total_dist = total_dist + seg_dist;
    %path_spd = (uv_map(i,j)+uv_map(i-1,j+1))/2; %point to point speed, avg of start and endpoint (kts)
    %path_spd = uv_map_lead(i,j); %simplified to get around index error
    seg_time = seg_dist/path_spd; %hours
    map_time = map_time + seg_time; %hours
    total_time = total_time + seg_time; %hours
    i = i-1;
    startlat = startlat-1;
    j = j+1;
    startlong = startlong+1;
elseif next_box == 3
    max_index = 3;
    long_path(m,k) = long_map(startlat,startlong+1); %long of desired gradient
    lat_path(m,k) = lat_map(startlat,startlong+1); %lat of desired gradient
    seg_dist = SegmentDistance(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(startlat,startlong+1),long_map(startlat,startlong+1));
    total_dist = total_dist + seg_dist;
    %path_spd = (uv_map(i,j)+uv_map(i,j+1))/2; %point to point speed, avg of start and endpoint (kts)
    %path_spd = uv_map_lead(i,j); %simplified to get around index error
    seg_time = seg_dist/path_spd; %hours
    map_time = map_time + seg_time; %hours
    total_time = total_time + seg_time; %hours
    i = i;
    startlat = startlat;
    j = j+1;
    startlong = startlong+1;
elseif next_box == 4
    max_index = 4;
    long_path(m,k) = long_map(startlat+1,startlong+1); %long of desired gradient
    lat_path(m,k) = lat_map(startlat+1,startlong+1); %lat of desired gradient
    seg_dist = SegmentDistance(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(startlat+1,startlong+1),long_map(startlat+1,startlong+1));
    total_dist = total_dist + seg_dist;
    %path_spd = (uv_map(i,j)+uv_map(i+1,j+1))/2; %point to point speed, avg of start and endpoint (kts)
    %path_spd = uv_map_lead(i,j); %simplified to get around index error
    seg_time = seg_dist/path_spd; %hours
    map_time = map_time + seg_time; %hours
    total_time = total_time + seg_time; %hours
    i = i+1;
    startlat = startlat+1;
    j = j+1;
    startlong = startlong+1;
elseif next_box == 0
    max_index = 0;
    long_path(m,k) = long_map(startlat-1,startlong-1); %long of desired gradient
    lat_path(m,k) = lat_map(startlat-1,startlong-1); %lat of desired gradient
    seg_dist = SegmentDistance(lat_map(startlat,startlong),long_map(startlat,startlong),lat_map(startlat-1,startlong-1),long_map(startlat-1,startlong-1));
    total_dist = total_dist + seg_dist;
    %path_spd = (uv_map(i,j)+uv_map(i+1,j+1))/2; %point to point speed, avg of start and endpoint (kts)
    %path_spd = uv_map_lead(i,j); %simplified to get around index error
    seg_time = seg_dist/path_spd; %hours
    map_time = map_time + seg_time; %hours
    total_time = total_time + seg_time; %hours
    i = i-1;
    startlat = startlat-1;
    j = j-1;
    startlong = startlong-1;
end

end



