function desc=spatialGrid(img, ndiv)
    img = double(img)./255;
    height = size(img, 1);
    width = size(img, 2);
    windows_h = ceil(height/ndiv);
    windows_w = ceil(width/ndiv);
    desc = [];
    for y=0:ndiv-1
        for x=0:ndiv-1
            offset_y = y*windows_h;
            offset_x = x*windows_w;
            img_cell = img(offset_y + 1:min(offset_y+windows_h, height), ...
                offset_x+1:min(offset_x+windows_w, width), :);
            avg_red = mean(img_cell(:, :, 1), 'all');
            avg_green = mean(img_cell(:, :, 2), 'all');
            avg_blue = mean(img_cell(:, :, 3), 'all');
            desc = [desc avg_red avg_green avg_blue];
        end
    end
end