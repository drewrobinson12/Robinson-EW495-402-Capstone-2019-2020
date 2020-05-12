function [v_w_map,u_w_map] = WindSizer(file,lat_map,long_map)
% Subset EDMAPS grid from 32S to 56.5N, 278.5 to 352
% This will produce a 64x199 matrix


u_w = ncread(file,'u-component_of_wind_height_above_ground');
v_w = ncread(file,'v-component_of_wind_height_above_ground');
lat_w = ncread(file, 'lat');
long_w = ncread(file, 'lon');
heights = ncread(file, 'height_above_ground1');
u_w = double(u_w(:,:,1));
v_w = double(v_w(:,:,1));


dims = size(u_w);
long_w2 = transpose(long_w);
long_w_map = [];
lat_w_map = [];
for i = 1:dims(2)
    long_w_map(i,:) = long_w2(1,:);
end 

for k = 1:dims(1)
    lat_w_map(:,k) = lat_w(:,1);
end

u_w_map_decimal = transpose(u_w);
u_w_map = u_w_map_decimal;

v_w_map_decimal = transpose(v_w);
v_w_map = v_w_map_decimal;
% updated size matcher
lat_w2 = lat_w_map(2:44,12:144);
long_w2 = long_w_map(2:44,12:144);
u_w2 = u_w_map(2:44,12:144);
v_w2 = v_w_map(2:44,12:144);

sizey = size(lat_map);

f=0;
for d = 3:3:sizey(1)
    lat_w3(d-2,:) = lat_w2(d-2-f,:);
    lat_w3(d-1,:) = lat_w2(d-1-f,:);
    lat_w3(d,:) = lat_w2(d-1-f,:);
    lat_w3(d+1,:) = lat_w2(d-f,:);
    f = f+1;
end

f=0;
for d = 3:3:sizey(2)
    long_w3(:,d-2) = long_w2(:,d-2-f);
    long_w3(:,d-1) = long_w2(:,d-1-f);
    long_w3(:,d) = long_w2(:,d-1-f);
    long_w3(:,d+1) = long_w2(:,d-f);
    f = f+1;
end

f=0;
for d = 3:3:sizey(1)
    u_w3(d-2,:) = u_w2(d-2-f,:);
    u_w3(d-1,:) = u_w2(d-1-f,:);
    u_w3(d,:) = u_w2(d-1-f,:);
    u_w3(d+1,:) = u_w2(d-f,:);
    v_w3(d-2,:) = v_w2(d-2-f,:);
    v_w3(d-1,:) = v_w2(d-1-f,:);
    v_w3(d,:) = v_w2(d-1-f,:);
    v_w3(d+1,:) = v_w2(d-f,:);
    f = f+1;
end

f=0;
for d = 3:3:sizey(2)
    u_w4(:,d-2) = u_w3(:,d-2-f);
    u_w4(:,d-1) = u_w3(:,d-1-f);
    u_w4(:,d) = u_w3(:,d-1-f);
    u_w4(:,d+1) = u_w3(:,d-f);
    v_w4(:,d-2) = v_w3(:,d-2-f);
    v_w4(:,d-1) = v_w3(:,d-1-f);
    v_w4(:,d) = v_w3(:,d-1-f);
    v_w4(:,d+1) = v_w3(:,d-f);
    f = f+1;
end

for i = 1:sizey(1)
        long_w4(i,:) = long_w3(1,:);
    end
    
for k = 1:sizey(2)
        lat_w4(:,k) = lat_w3(:,1);
end

v_w_map = v_w4;
u_w_map = u_w4;

end

