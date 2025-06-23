function [vertices, Bmin, Bmax] = generate_grid_vertices(lon_min, lon_max, lat_min, lat_max, h_min, h_max)
    % 生成网格的8个顶点坐标
    lats = [lat_min, lat_max];
    lons = [lon_min, lon_max];
    hs = [h_min, h_max];
    
    vertices = [];
    for i = 1:2
        for j = 1:2
            for k = 1:2
                % [x, y, z] = blh2xyz(lats(i), lons(j), hs(k));
                [x, y, z] = BLHtoXYZ_sphere2(lats(i), lons(j), hs(k));
                vertices = [vertices; x y z];
            end
        end
    end
    
    % 计算AABB的Bmin和Bmax
    Bmin = min(vertices);
    Bmax = max(vertices);
end