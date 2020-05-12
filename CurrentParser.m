function [lat_map,long_map,u_map,v_map,rows,cols] = CurrentParser(map_set)
% Subsetting 35S to 56N and 284 to 350 makes a 64x199 matrix
% Parses out lat, long, u, v of current and matrix dims

    for q = 1:length(map_set)
        file = 'oscar_vel'+map_set(q)+'.nc_subset.nc';
        lat_nc = ncread(file,'latitude');
        long_nc = ncread(file,'longitude');
        u_nc = ncread(file,'u');
        v_nc = ncread(file,'v');

        dims = size(u_nc);
        long_nc2 = transpose(long_nc);

        for i = 1:dims(2)
            long_map(i,:) = long_nc2(1,:);
        end

        for k = 1:dims(1)
            lat_map(:,k) = lat_nc(:,1);
        end

        u_map(:,:,q) = transpose(u_nc);
        v_map(:,:,q) = transpose(v_nc);
        %uv_map(:,:,q) = hypot(u_map(:,:,q),v_map(:,:,q));
        [rows, cols, depth] = size(u_map);

    end
end

