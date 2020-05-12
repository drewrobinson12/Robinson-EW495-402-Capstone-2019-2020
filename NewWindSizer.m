%% New Wind Sizer


file = 'uwind_19.nc'

u_test = ncread(file,'uwnd');
%v_test = ncread(file,);
lat_test = ncread(file, 'lat');
long_test = ncread(file, 'lon');

u_test_1 = transpose(u_test(:,:,1))
long_test1 = long_test(3:end);
lat_test1 = lat_test(1:end);


%% 

lat_hi = lat_test(1);
lat_lo = lat_test(end-1);

long_lo = long_test(3);
long_hi = long_test(end);

lat_test2 = [];
for i = 2:12
    lat_test2(i,1) = (lat_test1(i-1,1)+lat_test1(i+1,1))/2
end

%% New try interp2
%matched to: https://www.mathworks.com/matlabcentral/answers/113230-how-to-regrid-interpolate-netcdf-data-2-5x2-5-degree-to-higher-resolution-e-g-0-25x0-25-degree

[OY,OX] = meshgrid(283.125:1.875:350.625, 35.2375:1.9047:56.1893) %this is 12x37
[NY,NX] = meshgrid(283.125:0.34:350.625, 35.2375:0.33:56.1893) %This makes 64x199

u_test2 = u_test_1(1:12,1:37);
new_u = interp2(OY,OX,u_test2,NY,NX)

%% 2019 u-data
%matched to: https://www.mathworks.com/matlabcentral/answers/113230-how-to-regrid-interpolate-netcdf-data-2-5x2-5-degree-to-higher-resolution-e-g-0-25x0-25-degree

[OY,OX] = meshgrid(283.125:1.875:350.625, 35.2375:1.9047:56.1893) %this is 12x37
[NY,NX] = meshgrid(283.125:0.34:350.625, 35.2375:0.33:56.1893) %This makes 64x199

for i = 1:609
    u_trim = transpose(u_test(1:37,1:12,i)); %raw davies data, trimmed for correct bounds and transposed. 39x13 to 12x37.
    uw_19(:,:,i) = interp2(OY,OX,u_trim,NY,NX);
end

%% 2019 v-data

file = 'vwind_19.nc'

v_raw19 = ncread(file,'vwnd');

for i = 1:609
    v_trim = transpose(v_raw19(1:37,1:12,i)); %raw davies data, trimmed for correct bounds and transposed. 39x13 to 12x37.
    vw_19(:,:,i) = interp2(OY,OX,v_trim,NY,NX);
end

%% 2018 u-data
file = 'uwind_18.nc'

u_raw18 = ncread(file,'uwnd');

for i = 1:609
    u_trim = transpose(u_raw18(1:37,1:12,i)); %raw davies data, trimmed for correct bounds and transposed. 39x13 to 12x37.
    uw_18(:,:,i) = interp2(OY,OX,u_trim,NY,NX);
end

%% 2018 v-data

file = 'vwind_18.nc'

v_raw18 = ncread(file,'vwnd');

for i = 1:609
    v_trim = transpose(v_raw18(1:37,1:12,i)); %raw davies data, trimmed for correct bounds and transposed. 39x13 to 12x37.
    vw_18(:,:,i) = interp2(OY,OX,v_trim,NY,NX);
end
