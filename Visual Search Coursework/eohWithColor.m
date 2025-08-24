function desc=eohWithColor(img, ndiv, q_level, threshold)
    norm_img = double(img)./255;
    grey_img = norm_img(:, :, 1)*0.30 + norm_img(:, :, 2)*0.59 + norm_img(:, :, 3) * 0.11;
    grey_img = imgaussfilt(grey_img);
    filter_x = [1 0 -1; 2 0 -2; 1 0  -1]./4;
    filter_y = filter_x.';
    height = size(grey_img, 1);
    width = size(grey_img, 2);
    windows_h = ceil(height/ndiv);
    windows_w = ceil(width/ndiv);
    desc = [];
    for y=0:ndiv-1
        for x=0:ndiv-1
            offset_y = y*windows_h;
            offset_x = x*windows_w;
            img_cell = grey_img(offset_y + 1:min(offset_y+windows_h, height), ...
                offset_x+1:min(offset_x+windows_w, width));
            dx = conv2(img_cell, filter_x, 'same'); 
            dy = conv2(img_cell, filter_y, 'same');
            theta = rad2deg(atan2(dy, dx) + pi)/361; %[-pi, pi]->[0, 2pi]->[0, 360]->[0, 1)
            mag = sqrt (dx.^2 + dy.^2);
            filtered = theta(mag > threshold);
            bin = floor(filtered .* q_level);
            bin = reshape(bin, 1, []);
            tex_desc = histogram(bin, q_level).Values;
            total_freq = sum(tex_desc);
            if total_freq
                tex_desc = tex_desc ./ sum(tex_desc);
            end

            img_cell = norm_img(offset_y + 1:min(offset_y+windows_h, height), ...
                offset_x+1:min(offset_x+windows_w, width), :);
            avg_red = mean(img_cell(:, :, 1), 'all');
            avg_green = mean(img_cell(:, :, 2), 'all');
            avg_blue = mean(img_cell(:, :, 3), 'all');
            desc = [desc tex_desc avg_red avg_green avg_blue];
        end
    end
end