%% Multitrack iterative lead/lag

%% EDMAPS
%days since October 12, 1992
%top row is 24May18 to 24Aug18, bottom is 14May19 to 08Sep19
maps2018 = string([9355 9360 9366 9371 9376 9381 9386 9391 9396 9401 9406 9411 9416 9421 9426 9431 9436 9442 9447]);
maps2019 = string([9710 9715 9720 9725 9731 9736 9741 9746 9751 9756 9761 9766 9771 9776 9781 9786 9791 9796 9801 9807 9812 9817 9822 9827]);
 
%% Initial Clear
long_map = [];
lat_map = [];
u_map = [];
v_map = [];
uv_map_lead = [];

%Start and finish microtransat
load coastlines
startx = [313 313 295 283];
starty = [48 45.5 40 30];
finishx = [344 344 349 346.5 346.5];
finishy = [55 51 46 44.5 37];

%% Begin for loop of collecting current lat, long, u, v from each dataset
%ensure map set is for correct year
% Subsetting 35S to 56N and 284 to 350 makes a 64x199 matrix
% Parses out lat, long, u, v of current and matrix dims

[lat_map,long_map,u_map,v_map,rows,cols] = CurrentParser(maps2018); %change row 236 as well
    
%% Matching wind data size to current size
% Subset EDMAPS grid from 32S to 56.5N, 278.5 to 352
% This will produce a 64x199 matrix
%file = '19Aug19_wind.nc';
%[v_w_map,u_w_map] = WindSizer(file,lat_map,long_map);

%Upload wind sets from storage
load('WindData_FINAL_14APR.mat')

%to keep convention, BUT CHECK YEAR
v_w_map = vw_18;
u_w_map = uw_18;

%if help/understanding needed, look through NewWindSizer function

%% Post-Processing --Statistics
for i = 1:rows
    for j = 1:cols
        means_u(i,j) = mean(u_map(i,j,:));
        means_v(i,j) = mean(v_map(i,j,:));
        stds_u(i,j) = std(u_map(i,j,:));
        stds_v(i,j) = std(v_map(i,j,:));
    end
end 
%% Everything from all the above sections in one WS
load('InitWorkspace.mat')

%delete coastlat and long, and means and stds

%% 2018 Trimmed
load('InitWS2018.mat')
%% Everything except big unnecessary stuff

load('InitWorkspaceTrimmed.mat') %2019
%% Iterative segment mapper

%In this, I have altered the standard 1:rows/cols with variables that
%represent segmentation
lat_range = [];
long_range = [];

long_lookout_w = 1; % # gridsquares to project out to per segment TO THE EAST, every 3 = 1 degree
long_lookout_e = 9;
startlong = 3; %IMPORTANT** INITIAL START COLUMN
startlong_lag = 3;
endlong_w = startlong - long_lookout_w;
endlong_e = startlong + long_lookout_e - 1;
long_lookout = long_lookout_w + long_lookout_e;

lat_lookout_n = 6; % # gridsquares to project out to per segment TO THE NORTH, every 3 = 1 degree,
lat_lookout_s = 1; % # gridsquares to project out to per segment TO THE SOUTH, every 3 = 1 degree,
startlat = 63; %IMPORTANT** INITIAL START ROW (61 lowest w/out index error)
startlat_lag = 63;
endlat_n = startlat - lat_lookout_n + 1;
endlat_s = startlat + lat_lookout_s;
lat_diff = endlat_s-endlat_n;
lat_lookout = lat_lookout_n + lat_lookout_s;

lat_range = endlat_n:endlat_s;
long_range = endlong_w:endlong_e;

n_map = 2; %index of current map to use
%n_wind initialized in p = 1:r
n = 1; %number of trials
k = 1; %while iteration counter
c = 1; %plotting row counter
ct = 0; %counter for how many times wind maps recycle
min_count = 0; %count under current threshold
k1 = 1;

b2t_lat = [16,22,31,36,47,58];
b2t_long = [181,189,196,189,189];
%tgt_idx = 2;
total_time = 0;
total_dist = 0;
total_time_lag = 0;
total_dist_lag = 0;
map_time = 0;
map_time_lag = 0;

%sims = [];
long_path = [];
lat_path = [];
lat_path_lag = [];
long_path_lag = [];
%long_path_lag_plot = zeros(20);

for tgt_idx = 1
    success_cnt = 0;
    success_rate = 0;
    tgt_dist = 0;
    %clear sim %clear previous stuff
    %clear '29Apr20_Tgt1_5000_3.mat' %clear previous stuff
    name = 'test4';
    r = 2; %number of simulations
    inc = 1;
    iter = 1;
    %tgt_idx = 1;
    flag = 0;
    finish_fail = 0;
    twocount = 0;
    twocountlag = 0;
    onecount = 0;
    onecountlag = 0;
    agg_time = 0;
    agg_dist = 0;
    agg_spd = 0;
    agg_wind_maps = 0;
    agg_current_maps = 0;
%     agg_dist_array = [];
%     agg_time_array = [];
%     agg_spd_array = [];
%     agg_wind_maps_array = [];
%     agg_current_maps_array = [];
%     agg_tgt_dist_array = [];


    for p = 1:r
        
        lat_range = [];
        long_range = [];
        
        long_lookout_w = 0; % # gridsquares to project out to per segment TO THE EAST, every 3 = 1 degree
        long_lookout_e = 9;
        startlong = 4; %IMPORTANT** INITIAL START COLUMN
        startlong_lag = 4;
        endlong_w = startlong - long_lookout_w;
        endlong_e = startlong + long_lookout_e - 1;
        long_lookout = long_lookout_w + long_lookout_e;
        
        lat_lookout_n = 6; % # gridsquares to project out to per segment TO THE NORTH, every 3 = 1 degree,
        lat_lookout_s = 1; % # gridsquares to project out to per segment TO THE SOUTH, every 3 = 1 degree,
        startlat = 63; %IMPORTANT** INITIAL START ROW (61 lowest w/out index error)
        startlat_lag = 63;
        endlat_n = startlat - lat_lookout_n + 1;
        endlat_s = startlat + lat_lookout_s;
        lat_diff = endlat_s-endlat_n;
        lat_lookout = lat_lookout_n + lat_lookout_s;
        
        lat_range = endlat_n:endlat_s;
        long_range = endlong_w:endlong_e;
        
        n_map = 2; %index of current map to use
        n_wind = randi([3, 609]); %index of wind map to use
        %n_wind = 5;
        n = 1; %number of trials
        flag = 0;
        wind_cycles = 0;
        current_cycles = 0;
        k = 1; %while iteration counter
        klag = 1;
        c = 1; %plotting row counter
        ct = 0; %counter for how many times wind maps recycle
        min_count = 0; %count under current threshold
        k1 = 1;
        dist = 0;
        
        b2t_lat = [16,22,31,36,47,58];
        b2t_long = [181,189,196,189,189];
        total_time = 0;
        total_dist = 0;
        total_time_lag = 0;
        total_dist_lag = 0;
        map_time = 0;
        map_time_lag = 0;
        
        sims = [];
        long_path = [];
        lat_path = [];
        lat_path_lag = [];
        long_path_lag = [];
        %long_path_lag_plot = zeros(20);
        
        
        
        %figure(p);clf
        %figure(400);clf
        %fig = figure('Visible','off'); clf
        datapoint_u = [];
        datapoint_v = [];
        data_u_lag = [];
        data_v_lag = [];
        generated_uv_lead = [];
        generated_uv_lag = [];
        curr_angle = [];
        curr_angle_lag = [];
        
        finish_interx = [];
        finish_intery = [];
        start_interx = [];
        start_intery = [];
        %begin simulation for loop
        for m = 1:n %eventually up to 1000
            while lat_range(1) >= 1 && lat_range(end) <= 64 && long_range(end) < cols && n_map <= length(maps2018) %keeps within bounds of global N/S/E/W boundaries and 100 days
                if flag == 1
                    onecount = onecount+1;
                    total_dist = 0;
                    total_time = 0;
                    break
                elseif flag == 2
                    p = p-1;
                    twocount = twocount+1;
                    total_dist = 0;
                    total_time = 0;
                    break
                end
                %         %probably need to revisit to get this to update more often
                %         wind_counter = floor(map_time/6);
                %         n_wind = n_wind + wind_counter;
                %         wind_cycles = wind_cycles + 1;
                %         if n_wind > 609
                %             ct = ct+1;
                %             n_wind = randi([5, 609]);
                %         end
                %         n_wind_lag = n_wind - 4; %this will put it a full day behind
                
                if map_time > 70
                    map_time = 0;
                    if n_map == length(maps2018)
                        n_map = 2;
                        current_cycles = current_cycles + 1;
                    else
                        n_map = n_map + 1;
                        current_cycles = current_cycles + 1;
                    end
                    %klag = klag+1;
                    %lat_path_lag(1,klag) = NaN;
                    %long_path_lag(1,klag) = NaN;
                    %klag = klag+1;
                end
                
                %keep the wind map tracking
                
                
                
                
                
                endlong_w = startlong - long_lookout_w;
                endlong_e = startlong + long_lookout_e - 1;
                long_lookout = long_lookout_w + long_lookout_e;
                
                endlat_n = startlat - lat_lookout_n + 1;
                endlat_s = startlat + lat_lookout_s;
                lat_diff = endlat_s-endlat_n;
                lat_lookout = lat_lookout_n + lat_lookout_s;
                
                
                
                if endlat_s > 64
                    endlat_s = 64; %trying to prevent out of bounds index errors (64 is global row max)
                end
                if endlong_e > 199
                    endlong_e = 199;
                end
                if endlat_n < 1
                    endlat_n = 1;
                end
                
                lat_range = endlat_n:endlat_s;
                long_range = endlong_w:endlong_e;
                
                for i = 1:length(lat_range)
                    ci = lat_range(i);
                    for j = 1:length(long_range)
                        cj = long_range(j);
                        %datapoint_u(i,j,m) = 1.94*random('Normal',means_u(ci,cj),stds_u(ci,cj)); %insert these when you want "random" generate maps
                        %datapoint_v(i,j,m) = 1.94*random('Normal',means_v(ci,cj),stds_v(ci,cj));
                        datapoint_u(i,j,m) = 1.94*u_map(ci,cj,n_map); %insert these when you want true EDMAPS data maps
                        datapoint_v(i,j,m) = 1.94*v_map(ci,cj,n_map);
                        data_u_lag(i,j,m) = 1.94*u_map(ci,cj,n_map-1);
                        data_v_lag(i,j,m) = 1.94*v_map(ci,cj,n_map-1);
                        %                 datapoint_uw(i,j,m) = 1.94*u_w_map(ci,cj,n_wind); %wind, kts
                        %                 datapoint_vw(i,j,m) = 1.94*v_w_map(ci,cj,n_wind);
                        %                 data_uw_lag(i,j,m) = 1.94*u_w_map(ci,cj,n_wind_lag); %wind, kts
                        %                 data_vw_lag(i,j,m) = 1.94*v_w_map(ci,cj,n_wind_lag);
                        %Note: 1.94 is scaling factor converting m/s to kts
                        curr_angle(i,j,m) = unit2comp(atan2(datapoint_v(i,j,m),datapoint_u(i,j,m)));
                        curr_angle_lag(i,j,m) = unit2comp(atan2(data_v_lag(i,j,m),data_u_lag(i,j,m)));
                    end
                end
                
                generated_uv_lead = hypot(datapoint_u(:,:,m),datapoint_v(:,:,m));
                generated_uv_lag = hypot(data_u_lag(:,:,m),data_v_lag(:,:,m));
                %generated_uvw = hypot(datapoint_uw(:,:,m),datapoint_vw(:,:,m));
                %generated_uvw_lag = hypot(data_uw_lag(:,:,m),data_vw_lag(:,:,m));
                % ***Matrix Mapper**
                
                %Universal Declarations
                uv_map_lead = generated_uv_lead;
                uv_map_lag = generated_uv_lag;
%                 uv_w_map = generated_uvw;
%                 uv_w_map_lag = generated_uvw_lag;
                
                i = round(lat_lookout/2);
                j = 2;
                ii = round(lat_lookout/2);
                jj = 2;
                %k = 1;
                
                seg_dist = 0;
                seg_time = 0;
                path_spd = 0;
                
                check_area_init = [];
                check_area = [];
                curr_check_area_init = [];
                curr_check_area = [];
                heading_path = [];
                %long_path = [];
                %lat_path = [];
                %long_path(1) = long_map(1); %keeps from global zero plotting due to NaN
                %lat_path(1) = lat_map(end);
                %         startlong_hold = startlong;
                %         startlat_hold = startlat;
                %         startlong_lag = startlong_hold;
                %         startlat_lag = startlat_hold;
                %
                startlong_lag = startlong;
                startlat_lag = startlat;
                
                while i ~= find(lat_range==endlat_s) && i~= find(lat_range==endlat_n) && j ~= find(long_range==endlong_w) && j~= find(long_range==endlong_e) && ii ~= find(lat_range==endlat_s) && ii~= find(lat_range==endlat_n) && jj ~= find(long_range==endlong_w) && jj~= find(long_range==endlong_e) && seg_time <= 70
                    %Experimental wind update
             
             
                    %probably need to revisit to get this to update more often
                    wind_counter = floor(map_time/6);
                    n_wind = n_wind + wind_counter;
                    wind_cycles = wind_cycles + 1;
                    if n_wind > 609
                        ct = ct+1;
                        n_wind = randi([3, 609]);
                        %n_wind = 3;
                    end
                    n_wind_lag = n_wind - 2; %this will put it a half day behind
                    
                    for s = 1:length(lat_range)
                        ci = lat_range(s);
                        for t = 1:length(long_range)
                            cj = long_range(t);
                            
                            datapoint_uw(s,t,m) = 1.94*u_w_map(ci,cj,n_wind); %wind, kts
                            datapoint_vw(s,t,m) = 1.94*v_w_map(ci,cj,n_wind);
                            data_uw_lag(s,t,m) = 1.94*u_w_map(ci,cj,n_wind_lag); %wind, kts
                            data_vw_lag(s,t,m) = 1.94*v_w_map(ci,cj,n_wind_lag);
                        end
                    end
                    
                    k = k+1; %iteration counter
                    klag = klag+1;
                    n_wpts = k-1; %waypoint counter
                    n_segs = k-2; %segment counter
                    
                    result_u = [];
                    result_v = [];
                    result_mag = [];
                    result_dir = [];
                    wsidot_u = [];
                    wsidot_v = [];
                    curr_u = [];
                    curr_v = [];
                    
                    %**** START OF LEAD IF*********
                    %CURRENT ONLY
                    %[lat_path,long_path,seg_dist,total_dist,path_spd,seg_time,total_time,i,j,startlat,startlong] = LeadIf(k,m,i,j,lat_range,endlat_s,endlat_n,long_range,endlong_e,endlong_w,uv_map_lead,curr_angle,startlat,startlong,lat_map,long_map,total_dist,total_time,lat_path,long_path,b2t_lat,b2t_long,tgt_idx);
                    %         startlong_hold = startlong;
                    %         startlat_hold = startlat;
                    %COMBINED WIND AND CURRENT
                    [lat_path,long_path,seg_dist,total_dist,path_spd,seg_time,total_time,i,j,startlat,startlong,min_count,map_time,flag] = LeadIf2(k,m,i,j,lat_range,endlat_s,endlat_n,long_range,endlong_e,endlong_w,uv_map_lead,curr_angle,startlat,startlong,lat_map,long_map,total_dist,total_time,lat_path,long_path,b2t_lat,b2t_long,tgt_idx,datapoint_uw,datapoint_vw,datapoint_u,datapoint_v,min_count,seg_dist,seg_time,map_time);
                    if flag == 1 %boundary error, want to count as failure
                        %onecount = onecount+1;
                        total_dist = 0;
                        total_time = 0;
                        break
                    elseif flag == 2 %initial corner area, want to redo/not count towards success or failure
                        %twocount = twocount+1;
                        total_dist = 0;
                        total_time = 0;
                        break
                    end
                    
                    
                    %**** START OF LAG IF********
                    %         startlong_lag = startlong_hold;
                    %         startlat_lag = startlat_hold;
                    %long_path_lag = [];
                    %lat_path_lag = [];
                    [lat_path_lag,long_path_lag,seg_dist,total_dist_lag,path_spd,seg_time,total_time_lag,ii,jj,startlat_lag,startlong_lag,min_count,map_time_lag,flag] = LagIf2(klag,m,ii,jj,lat_range,endlat_s,endlat_n,long_range,endlong_e,endlong_w,uv_map_lag,curr_angle_lag,startlat_lag,startlong_lag,lat_map,long_map,total_dist_lag,total_time_lag,lat_path_lag,long_path_lag,b2t_lat,b2t_long,tgt_idx,data_uw_lag,data_vw_lag,datapoint_u,datapoint_v,min_count,seg_dist,seg_time,map_time_lag);
                    if flag == 1 %boundary error
                        onecountlag = onecountlag+1;
                        total_dist_lag = 0;
                        total_time_lag = 0;
                        break
                    elseif flag == 2 %initial corner area, want to redo/not count towards success or failure
                        twocountlag = twocountlag + 1;
                        break
                    end
                    
                end %end local while loop
                %{
        if map_time > 120
            lat_path_lag(1,end) = NaN;
            long_path_lag(1,end) = NaN;
            k = k+1;
        end
                %}
                
                
                %{
       if mean(long_map(1,startlong)) > 298
            break;
        end
                %}
                
                
                %**Plotting full view** DO NOT use both full view and segment
                %             for b = 1:c
                %             long_path_lag_plot(b,:) = long_path_lag(b,k1:k);
                %             end
                %{
            figure(400);
            
            %subplot(2,2,1)
%                grid on
%                plot(long_path(long_path~=0),lat_path(lat_path~=0)) %remove zeros
%                hold on
%                plot(long_path_lag(long_path_lag~=0),lat_path_lag(lat_path_lag~=0))
%                hold on
                
            %

       %figure(p); %regular plot
       %grid on
       %{ %this is for plotting MicroT lines
       %
            %start line
            plot(startx,starty,'g','LineWidth',2)
            hold on
            %finish line
            plot(finishx,finishy,'r','LineWidth',2)
            hold on
         %
       %
            title('Actual/Uncertain/Random Route')
        
            xlabel("Longitude (" + char(176) + "E)")
            ylabel("Latitude (" + char(176) + "N)")
        
            plot((coastlon+360),coastlat,'b')
            daspect([2 2 2])
            grid on
            %axesm('MapProjection','mercator','MapLonLimit',[-80 0],'MapLatLimit',[30 60])
            %axis([282 350 33 50])
            %axis([long_map(startlong)-2 long_map(1,endlong)+2 lat_map(endlat+lat_diff,1)-2 lat_map(endlat)+2])

            hold on
            %legend('Coastline','Lead/Actual','Lag/Predicted')
           
            %figure(400); %quiver current plot
                %{
            %subplot(2,2,2)
            title('Actual/Uncertain/Random Current')
            starthold = endlong_e-endlong_w;
            %startlat = endlat-lat_diff;
            grid on
            quiverC2D(long_map(endlat_n:endlat_s,endlong_w:endlong_e),lat_map(endlat_n:endlat_s,endlong_w:endlong_e),datapoint_u(1:length(lat_range),1:length(long_range),1),datapoint_v(1:length(lat_range),1:length(long_range),1))
            hold on
            daspect([2 2 2])
            axis([284 350 34 56])
            %axesm mercator
            %
            subplot(2,2,3)
            title('Real Current')
            quiverC2D(long_map,lat_map,u_map(:,:,n_map),v_map(:,:,n_map),60000)
            grid on
            hold on
            daspect([2 2 2])
            axis([284 350 34 56])
            %axesm mercator
                %}
            %subplot(2,2,4)
            title('Overlay')
            quiverC2D(long_map,lat_map,1.94*u_map(:,:,n_map),1.94*v_map(:,:,n_map),60000)
            grid on
            hold on
            plot(long_path(long_path~=0),lat_path(lat_path~=0),'k') %remove zeros
            hold on
            %plot(long_path_lag,lat_path_lag,'r')
            %plot(long_path_lag(long_path_lag~=0),lat_path_lag(lat_path_lag~=0),'r')
            hold on
            daspect([2 2 2])
            axis([284 350 34 56])
            
                %}
                %endlat_s = 1;
                c = c+1;
                %avg_spd = total_dist/total_time;
            end %end map boundary while
            %plotting
            %
            
            %fig
        %
            figure(400)
            
            %figure('Color','c');
            plot(startx,starty,'Color',[0.24,0.72,0.22],'LineWidth',2)
            hold on
            %finish line
            plot(finishx,finishy,'r','LineWidth',2)
            hold on
            
            %title('Actual/Uncertain/Random Route')
            
            xlabel("Longitude (" + char(176) + "W)")
            ylabel("Latitude (" + char(176) + "N)")
            xticklabels({'80','70','60','50','40','30','20','10'})
            
            grid on
            %title('USNA Sailbot - Microtransat Challenge')
            hold on
            %quiverC2D(long_map,lat_map,1.94*u_map(:,:,n_map),1.94*v_map(:,:,n_map),60000)
            %quiver(long_map,lat_map,1.94*u_map(:,:,n_map-1),1.94*v_map(:,:,n_map-1))
            %grid on
            %hold on
            plot(long_path(long_path~=0),lat_path(lat_path~=0),'b','LineWidth',0.5) %plot lead
            hold on
            %plot(long_path_lag,lat_path_lag,'r')
            plot(long_path_lag(long_path_lag~=0),lat_path_lag(lat_path_lag~=0),'m') %plot lag
            hold on
            daspect([2 2 2])
            axis([283 352 34 56])
            set(gca,'Color',[0.67,0.92,0.89]) %ocean color
            
            legend('Start Line','Finish Line','Actual Route (Lead)','Projected Route (Lag)','NumColumns',1,'AutoUpdate','off','Location','best','FontSize',10)
            hold on
            %plot((coastlon+360),coastlat,'b')
            mapshow(coastlon+360, coastlat, 'DisplayType','polygon','FaceColor',[0.07,0.4,0]) %land color
            daspect([2 2 2])
            
            %}
            %
            %**Determining Success**
            start_interx = polyxpoly(startx,starty,long_path(m,2:end),lat_path(m,2:end));
            [finish_interx,finish_intery] = polyxpoly(finishx,finishy,long_path(m,2:end),lat_path(m,2:end));
            
            if isempty(start_interx) == 0 && isempty(finish_interx) == 0
                success_cnt = success_cnt + 1;
                dist = DistanceCalculator(finish_intery(1),finish_interx(1),lat_map(b2t_lat(tgt_idx),1),long_map(1,b2t_long(tgt_idx)));
                
                tgt_dist = tgt_dist + dist; %running total of error from tgt index
                %avg_tgt_dist = tgt_dist/success_cnt; %average error (nm) from tgt index
                
                agg_time = agg_time + total_time; %total time of the entire simulation
                %avg_time = (agg_time/24)/success_cnt; %avg journey time of successful voyages, in hours
                
                avg_spd = total_dist/total_time;
                agg_spd = agg_spd + avg_spd; %sum of each simulations average speed
                %avg_sim_spd = agg_spd/success_cnt;
                
                agg_dist = agg_dist + total_dist;
                %avg_dist = agg_dist/success_cnt; %avg total distance traveled in successful voyages, nm
                
                agg_wind_maps = agg_wind_maps + wind_cycles;
                %avg_wind_maps = agg_wind_maps/success_cnt;
                
                agg_current_maps = agg_current_maps + current_cycles;
                %avg_current_maps = agg_current_maps/success_cnt;
            elseif isempty(start_interx) == 1 || isempty(finish_interx) == 1
                finish_fail = finish_fail + 1;
            end
            
            %viablesims = inc-twocount; 
            %success_rate = success_cnt/(viablesims) * 100; %success rate out of non-corner error trials
            
            if mod(p,inc) == 0
                  
%                     sim = matfile(name,'Writable',true);
%                     sim.viablesims(2,loop) = viablesims;
%                     sim.success_cnt(2,loop) = success_cnt;
%                     sim.finish_fail(2,loop) = finish_fail;
%                     sim.onecount(2,loop) = onecount;
%                     sim.twocount(2,loop) = twocount;
%                     sim.agg_time(2,loop) = agg_time;
%                     sim.agg_spd(2,loop) = agg_spd;
%                     sim.agg_dist(2,loop) = agg_dist;
%                     sim.tgt_dist(2,loop) = tgt_dist;
%                     sim.agg_current_maps(2,loop) = agg_current_maps;
%                     sim.agg_wind_maps(2,loop) = agg_wind_maps;    
                  
                 
                
                loop = p/inc;
                sim = matfile(name,'Writable',true);
                sim.success_cnt(iter,loop) = success_cnt;
                sim.onecount(iter,loop) = onecount;
                sim.twocount(iter,loop) = twocount;
                sim.finish_fail(iter,loop) = finish_fail-onecount-twocount;
                sim.agg_time(iter,loop) = agg_time;
                sim.agg_spd(iter,loop) = agg_spd;
                sim.agg_dist(iter,loop) = agg_dist;
                sim.tgt_dist(iter,loop) = tgt_dist;
                sim.agg_current_maps(iter,loop) = agg_current_maps;
                sim.agg_wind_maps(iter,loop) = agg_wind_maps;
                sim.viablesims(iter,loop) = success_cnt + onecount + (finish_fail-onecount-twocount);
                
                viablesims = 0;
                success_cnt = 0;
                finish_fail = 0;
                onecount = 0;
                twocount = 0;
                agg_time = 0;
                agg_spd = 0;
                agg_dist = 0;
                tgt_dist = 0;
                agg_current_maps = 0;
                agg_wind_maps = 0;
            end
            
        end
    end %end simulation loop
end

%%
    success_cnt = sum(sim.success_cnt,'all');
    viablesims = sum(sim.viablesims,'all');
    success_rate = success_cnt/viablesims * 100;
    onecount = sum(sim.onecount,'all');
    twocount = sum(sim.twocount,'all');
    finish_fail = sum(sim.finish_fail,'all');

    tgt_dist = sum(sim.tgt_dist,'all');
    avg_tgt_dist = tgt_dist/success_cnt; %average error (nm) from tgt index
    
    agg_time = sum(sim.agg_time,'all');
    avg_time = (agg_time/24)/success_cnt; %avg journey time of successful voyages, in hours
    
    agg_spd = sum(sim.agg_spd,'all');
    avg_sim_spd = agg_spd/success_cnt;
    
    agg_dist = sum(sim.agg_dist,'all');
    avg_dist = agg_dist/success_cnt; %avg total distance traveled in successful voyages, nm

    agg_wind_maps = sum(sim.agg_wind_maps,'all');
    avg_wind_maps = agg_wind_maps/success_cnt;
    
    agg_current_maps = sum(sim.agg_current_maps,'all');
    avg_current_maps = agg_current_maps/success_cnt;
    
    disp("Results of " + num2str(p) + " Simulations at Target Index " + num2str(tgt_idx))
    disp("Number of Viable Sims: " + num2str(viablesims))
    disp("Number of Successful Voyages: " + num2str(success_cnt) + ", or " + num2str(success_rate) + "%")
    disp("Number of Finish Line Failures: " + num2str(finish_fail-onecount-twocount))
    disp("Number of Boundary Failures: " + num2str(onecount))
    disp("Number of Corner Failures: " + num2str(twocount))
    disp("Average Time of Successful Voyages (days): " + num2str(avg_time))
    disp("Average Speed of Successful Voyages (kts): " + num2str(avg_sim_spd))
    disp("Average Distance Traveled of Successful Voyages (NM): " + num2str(avg_dist))
    disp("Average Distance from Target Lat/Long Coordinate (NM): " + num2str(avg_tgt_dist))
    disp("Average Current Map Update (days): " + num2str(avg_time/avg_current_maps))
    disp("Average Wind Map Update (days): " + num2str(avg_time/avg_wind_maps))
%end
%%
%Ocean Mapping
n_map = 5;

figure(401)
subplot(2,1,1)
 %figure('Color','c');
            plot(startx,starty,'Color',[0.24,0.72,0.22],'LineWidth',2)
            hold on
            %finish line
            plot(finishx,finishy,'r','LineWidth',2)
            hold on

            title('Current View')

            xlabel("Longitude (" + char(176) + "W)")
            ylabel("Latitude (" + char(176) + "N)")
            xticklabels({'80','70','60','50','40','30','20','10'})
            
            grid on
            %title('USNA Sailbot - Microtransat Challenge')
            hold on
            %quiverC2D(long_map,lat_map,1.94*u_map(:,:,n_map),1.94*v_map(:,:,n_map),60000)
            quiverC2D(long_map,lat_map,1.94*u_map(:,:,n_map),1.94*v_map(:,:,n_map),12736)
            %grid on
            %hold on
            %plot(long_path(long_path~=0),lat_path(lat_path~=0),'b','LineWidth',0.5) %plot lead
            hold on
            %plot(long_path_lag,lat_path_lag,'r')
            %plot(long_path_lag(long_path_lag~=0),lat_path_lag(lat_path_lag~=0),'m') %plot lag
            hold on
            daspect([2 2 2])
            axis([283 352 34 56])
            set(gca,'Color',[0.67,0.92,0.89]) %ocean color
            
            %legend('Start Line','Finish Line','Actual Route (Lead)','Projected Route (Lag)','NumColumns',1,'AutoUpdate','off','Location','eastoutside','FontSize',10)
            hold on
            %plot((coastlon+360),coastlat,'b')
            mapshow(coastlon+360, coastlat, 'DisplayType','polygon','FaceColor',[0.07,0.4,0]) %land color
            daspect([2 2 2])
            subplot(2,1,2)
            %figure('Color','c');
            plot(startx,starty,'Color',[0.24,0.72,0.22],'LineWidth',2)
            hold on
            %finish line
            plot(finishx,finishy,'r','LineWidth',2)
            hold on

            title('Current View')

            xlabel("Longitude (" + char(176) + "W)")
            ylabel("Latitude (" + char(176) + "N)")
            xticklabels({'80','70','60','50','40','30','20','10'})
            
            grid on
            %title('USNA Sailbot - Microtransat Challenge')
            hold on
            %quiverC2D(long_map,lat_map,1.94*u_map(:,:,n_map),1.94*v_map(:,:,n_map),60000)
            quiver(long_map,lat_map,1.94*u_map(:,:,n_map),1.94*v_map(:,:,n_map),'b')
            %grid on
            %hold on
            %plot(long_path(long_path~=0),lat_path(lat_path~=0),'b','LineWidth',0.5) %plot lead
            hold on
            %plot(long_path_lag,lat_path_lag,'r')
            %plot(long_path_lag(long_path_lag~=0),lat_path_lag(lat_path_lag~=0),'m') %plot lag
            hold on
            daspect([2 2 2])
            axis([283 352 34 56])
            set(gca,'Color',[0.67,0.92,0.89]) %ocean color
            
            %legend('Start Line','Finish Line','Actual Route (Lead)','Projected Route (Lag)','NumColumns',1,'AutoUpdate','off','Location','eastoutside','FontSize',10)
            hold on
            %plot((coastlon+360),coastlat,'b')
            mapshow(coastlon+360, coastlat, 'DisplayType','polygon','FaceColor',[0.07,0.4,0]) %land color
            daspect([2 2 2])


%% WIND PLOTTING - VISUAL AID
wind_map = 400;

figure(402);clf
title('Wind Map')
subplot(2,1,1)
plot(startx,starty,'Color',[0.24,0.72,0.22],'LineWidth',2)
hold on
%finish line
plot(finishx,finishy,'r','LineWidth',2)
hold on
quiverC2D(long_map,lat_map,100*u_w_map(:,:,wind_map),100*v_w_map(:,:,wind_map),12736)
grid on
hold on
daspect([2 2 2])
axis([283 352 34 56])
set(gca,'Color',[0.67,0.92,0.89]) %ocean color

%legend('Start Line','Finish Line','Actual Route (Lead)','Projected Route (Lag)','NumColumns',1,'AutoUpdate','off','Location','eastoutside','FontSize',10)
hold on
%plot((coastlon+360),coastlat,'b')
mapshow(coastlon+360, coastlat, 'DisplayType','polygon','FaceColor',[0.07,0.4,0]) %land color
daspect([2 2 2])

subplot(2,1,2)
title('Wind Map')
plot(startx,starty,'Color',[0.24,0.72,0.22],'LineWidth',2)
hold on
%finish line
plot(finishx,finishy,'r','LineWidth',2)
hold on
quiver(long_map,lat_map,1.94*u_w_map(:,:,wind_map),1.94*v_w_map(:,:,wind_map))
grid on
hold on
daspect([2 2 2])
axis([283 352 34 56])
set(gca,'Color',[0.67,0.92,0.89]) %ocean color

%legend('Start Line','Finish Line','Actual Route (Lead)','Projected Route (Lag)','NumColumns',1,'AutoUpdate','off','Location','eastoutside','FontSize',10)
hold on
%plot((coastlon+360),coastlat,'b')
mapshow(coastlon+360, coastlat, 'DisplayType','polygon','FaceColor',[0.07,0.4,0]) %land color
daspect([2 2 2])

%% Max wind speed
n_wind = 400;
max_spd = [];
for n_wind = 1:609
    max_spd(n_wind) = max(max(hypot(1.94*uw_19(:,:,n_wind),1.94*vw_19(:,:,n_wind))));
end

avg_max_speed = mean(max_spd)

mean_spd = [];
for n_wind = 1:609
    wspd = mean(hypot(1.94*uw_19(:,:,n_wind),1.94*vw_19(:,:,n_wind)),'all');
end

mean_wspd = mean(wspd)

%% Max curr speed
n_map = 1;
max_c_spd = [];
for n_map = 1:length(maps2019)
    max_c_spd(n_map) = max(max(hypot(1.94*u_map(:,:,n_map),1.94*v_map(:,:,n_map))));
end

avg_max_c_speed = mean(max_c_spd)

% mean_spd = [];
% for n_wind = 1:609
%     wspd = mean(hypot(uw_19(:,:,n_wind),vw_19(:,:,n_wind)),'all');
% end
% 
% mean_wspd = mean(wspd)


%% smoothing messing

lat_path_lag_s = fastsmooth(lat_path_lag,217,1,0);
long_path_lag_s = fastsmooth(long_path_lag,217,1,0);

figure(400);
            plot(startx,starty,'g','LineWidth',2)
            hold on
            %finish line
            plot(finishx,finishy,'r','LineWidth',2)
            hold on

            title('Actual/Uncertain/Random Route')

            xlabel("Longitude (" + char(176) + "E)")
            ylabel("Latitude (" + char(176) + "N)")

            
            grid on
            title('Overlay')
            hold on
            %quiverC2D(long_map,lat_map,1.94*u_map(:,:,n_map),1.94*v_map(:,:,n_map),60000)
            %quiver(long_map,lat_map,1.94*u_map(:,:,n_map-1),1.94*v_map(:,:,n_map-1))
            %grid on
            %hold on
            plot(long_path(long_path~=0),lat_path(lat_path~=0),'k','LineWidth',1) %remove zeros
            hold on
            %plot(long_path_lag,lat_path_lag,'r')
            plot(long_path_lag_s(long_path_lag_s~=0),lat_path_lag_s(lat_path_lag_s~=0),'m')
            hold on
            daspect([2 2 2])
            axis([284 350 34 56])
            legend('Start Line','Finish Line','Actual (Lead)','Projected (Lag)','NumColumns',1,'AutoUpdate','off','Location','eastoutside','FontSize',10)
            hold on
            plot((coastlon+360),coastlat,'b')
            daspect([2 2 2])

%% Distance Work
tot_dist = 0;
ivl = 6;
for i = 2:ivl:length(lat_path)-ivl
adj_dist = DistanceCalculator(lat_path(i),long_path(i),lat_path(i+ivl),long_path(i+ivl));
tot_dist = tot_dist + adj_dist;
end
format long g
disp(tot_dist)

%% Wind/Current intersect

curr_spd = [0:0.05:4];
wind_spd = [0:0.05:20];

boat_spd_w = 0.1175*wind_spd; %0.1175 = EV of scale factor coefficient
boat_spd_c = curr_spd;

figure(405);clf
title('Diff Axes')
yyaxis left
ylabel('Boat Speed from Current (kts)')
plot(curr_spd,boat_spd_c)

grid on
hold on
yyaxis right
ylabel('Boat Speed from Wind (kts)')
plot(wind_spd, boat_spd_w)

xlabel('Force Speed (kts)')

%% Ocean Colorbar

n_map = 4;

figure(407);clf

ax = gca;
ax.FontSize = 11
lab = xlabel("Longitude (" + char(176) + "W)")
lab.FontSize = 12
ylabel("Latitude (" + char(176) + "N)")

%title('USNA Sailbot - Microtransat Challenge')
hold on
colorbar('Ticks',[0 .25 0.5 0.75 1],'TickLabels',{'0 kts','1 kt','2 kts','3 kts','4 kts'},'FontSize',12)
quiverC2D(long_map,lat_map,u_map(:,:,n_map),v_map(:,:,n_map),12736)
%legend('Start Line','Finish Line','Actual Route (Lead)','Projected Route (Lag)','NumColumns',1,'AutoUpdate','off','Location','eastoutside','FontSize',10)
hold on
plot((coastlon+360),coastlat,'b')
%mapshow(coastlon+360, coastlat, 'DisplayType','polygon','FaceColor',[0.07,0.4,0]) %land color
axis([284 350 35 56])
xticklabels({'70','60','50','40','30','20','10'});
daspect([2 2 2])

%% Polar Chart Experimentation
labels = {90 60 30 0 330 300 270 240 210 180 150 120};
rh_labels = {0.10 0.11 0.12 0.13};

first_th = [0.8727:0.01:1.134]; %40 to 25
first_rho = 0.1*ones(size(first_th));

fst2_th = [2.007:0.01:2.269]; %335 to 320
fst2_rho = 0.1*ones(size(fst2_th));

third_th = [0.5236:0.01:0.8727]; %40 to 60
third_rho = 0.116*ones(size(third_th));

th2_th = [2.269:0.01:2.618]; %300 to 320
th2_rho = 0.116*ones(size(th2_th));

fourth_th = [6.109:0.01:6.4577]; %100 to 80
fourth_rho = 0.121*ones(size(fourth_th));

fo_splt_th = [6.4577:0.01:6.8068]; %80 to 60
fo_splt_rho = 0.123*ones(size(fo_splt_th));

fo2_th = [2.618:0.01:2.967]; %300 to 280
fo2_rho = 0.123*ones(size(fo2_th));

fo2_splt_th = [2.967:0.01:3.316]; %280 to 260
fo2_splt_rho = 0.121*ones(size(fo2_splt_th));

fifth_th = [5.411:0.01:6.109]; %140 to 100
fifth_rho = 0.116*ones(size(fifth_th));

fi2_th = [3.316:0.01:4.014]; %260 to 220
fi2_rho = 0.116*ones(size(fi2_th));

sixth_th = [5.0614:0.01:5.411];%160 to 140
sixth_rho = 0.119*ones(size(sixth_th));

si2_th = [4.014:0.01:4.363]; %220 to 200
si2_rho = 0.119*ones(size(si2_th));

second_th = [4.363:0.01:5.0614]; %200 to 160
second_rho = 0.123*ones(size(second_th));
% radial lines
seven_rho = [0.09:0.0001:0.1]%
seven_th = 1.134*ones(size(seven_rho))

eight_rho = [0.1:0.0001:0.116]%
eight_th = 0.8727*ones(size(eight_rho))

nine_rho = [0.116:0.0001:0.123]%[5:0.01:15]
nine_th = 6.8068*ones(size(nine_rho))

tween_rho = [0.121:0.0001:0.123];
tween_th = 6.4577*ones(size(tween_rho));

ten_rho = [0.116:0.0001:0.121]%[5:0.01:15]
ten_th = 6.109*ones(size(ten_rho))

elev_rho = [0.116:0.0001:0.119]%[5:0.01:15]
elev_th = 5.411*ones(size(elev_rho))

twel_rho = [0.119:0.0001:0.123]%[5:0.01:15]
twel_th = 5.0614*ones(size(twel_rho))

thtn_rho = [0.119:0.0001:0.123];
thtn_th = 4.363*ones(size(thtn_rho));

ftn_rho = [0.116:0.0001:0.119];
ftn_th = 4.014*ones(size(ftn_rho));

fftn_rho = [0.116:0.0001:0.121];
fftn_th = 3.316*ones(size(fftn_rho));

tween2_rho = [0.121:0.0001:0.123];
tween2_th = 2.967*ones(size(tween2_rho));

sxtn_rho = [0.116:0.0001:0.123];
sxtn_th = 2.618*ones(size(sxtn_rho));

svtn_rho = [0.1:0.0001:0.116];
svtn_th = 2.269*ones(size(svtn_rho));

eitn_rho = [0.09:0.0001:0.1];
eitn_th = 2.007*ones(size(eitn_rho));

figure(270);clf
ax = polaraxes;

polarscatter(first_th,first_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(fst2_th,fst2_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(second_th,second_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(third_th,third_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(th2_th,th2_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(fourth_th,fourth_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(fo_splt_th,fo_splt_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(fo2_th,fo2_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(fo2_splt_th,fo2_splt_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(fifth_th,fifth_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(fi2_th,fi2_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(sixth_th,sixth_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(si2_th,si2_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(seven_th,seven_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(eight_th,eight_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(nine_th,nine_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(tween_th,tween_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(ten_th,ten_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(elev_th,elev_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(twel_th,twel_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(thtn_th,thtn_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(ftn_th,ftn_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(fftn_th,fftn_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(tween2_th,tween2_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(sxtn_th,sxtn_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(svtn_th,svtn_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on
polarscatter(eitn_th,eitn_rho,'Marker','.','MarkerFaceColor','b','MarkerEdgeColor','b')
hold on

%ax = polaraxes
ax.RLim = [0.09 0.128];
rticks([0.10 0.11 0.12 0.13]);

thetaticklabels(labels);
rticklabels(rh_labels);%legend('Polar Curve')
%title("Experimental Polar Chart for Sailbot")
%rectangle('position',[0 0 x_nm y_nm])
%axis([0 x_nm 0 y_nm])

%% Success/failure plotting--Bar Graph

tgt = [1 2 3 4 5];
tgt_bkdown = [83.26917 14.55927 2.17155 ; 83.61914 13.80857 2.57228 ; 90.47478 6.284476 3.24074 ; 93.76418 2.97977 3.25604 ; 94.19075 1.47943 4.32981 ];
figure(500);clf
bar(tgt_bkdown,'stacked')
xlabel('Target Location')
ylabel('Rate (%)')
legend('Success','Finish Line Failure','Boundary Failure','Location','south')


%% Success/failure plotting--Pie Chart%
