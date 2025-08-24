function global_descriptors = siftDescriptors(image, num_intervals, num_octaves)
    image = imresize(image,4,"bicubic");
    norm_image = double(image)./255;
    grey_image = norm_image(:, :, 1)*0.30 + norm_image(:, :, 2)*0.59 + norm_image(:, :, 3) * 0.11;
    grey_image = imgaussfilt(grey_image, 1);
    grey_image = imsharpen(grey_image);
    image_height = size(grey_image, 1);
    image_width = size(grey_image, 2);
    if (num_intervals < 1)
        num_intervals = 2;
    end
    if (num_octaves < 1)
        num_octaves = floor(log10(min(image_height, image_width))/log10(2)) - 1;
    end
    k = 2 .^ (1/num_intervals);
    sigma_start = 1.6; 
    num_images_in_octave = num_intervals + 3;
    G = cell(num_octaves, 1); %Gaussian images
    sigmas = cell(num_octaves, 1);
    sigma = sigma_start;
    grey_image = padarray(grey_image, [1 1]);
    current_h = image_height+2;
    current_w = image_width+2;
    current_img = grey_image;
    for octave_ind=1:num_octaves
        octave_images = zeros(current_h, current_w, num_images_in_octave);
        current_sigmas = zeros(num_images_in_octave, 1);
        for level=1:num_images_in_octave
            octave_images(:, :, level) = imgaussfilt(current_img, sigma);
            current_sigmas(level) = sigma;
            sigma = k*sigma;
        end
        G{octave_ind, 1} = octave_images;
        sigmas{octave_ind, 1} = current_sigmas; %sigmas{octave_ind}(level)
        sigma = getfield(flip(current_sigmas), {3});
        current_img = current_img(2:2:end,2:2:end);
        current_h = floor(current_h/2);
        current_w = floor(current_w/2);
    end
    D ={}; %Difference of Gaussians
    for octave_ind=1:num_octaves
        diff_images = G{octave_ind}(:, :, 2:end) - G{octave_ind}(:, :, 1:end-1);
        D = [D; diff_images];
    end
    keypoints = [];
    for octave_ind=1:num_octaves
        for interval_ind=1:num_intervals
            extrema = local_extrema_detection(interval_ind, D{octave_ind}, octave_ind);
            keypoints = [keypoints; extrema];
        end
    end
    global_descriptors = aggregateDescriptors(G, keypoints, sigmas);
end

function global_descriptors=aggregateDescriptors(G, keypoints, sigmas)
    filter_x = [-1 1];
    filter_y = filter_x.';
    G_orientation = cell(size(G));
    G_magnitude = cell(size(G));
    sample_size = 16;
    desc_sigma = 1.5 * sample_size;
    sample_radius = 0.5 * sample_size;
    desc_array_size = 4;
    cell_size = sample_size/desc_array_size;
    exp_scale =  - 1/ (2*desc_sigma.^ 2);
    gauss_mag = exp_scale * (-1/pi);
    global_descriptors = [];
    for ind=1:size(keypoints, 1)
        % init
        octave_ind = keypoints(ind, 4);
        level = round(keypoints(ind, 3));
        sigma = 1.5 * sigmas{octave_ind}(level);

        %% create magnitude and orientation maps for gaussian image
        if isempty(G_orientation{octave_ind})
            gaussian_image = G{octave_ind}(:, :, level);
            gaussian_image = padarray(gaussian_image, [sample_radius sample_radius]);
            G_orientation{octave_ind} = zeros(size(gaussian_image, 1), size(gaussian_image, 2), size(G{octave_ind}, 3));
            G_magnitude{octave_ind} = zeros(size(G_orientation{octave_ind}));
            dx = conv2(gaussian_image, filter_x, 'same'); 
            dy = conv2(gaussian_image, filter_y, 'same');
            G_orientation{octave_ind}(:, :, level) = floor(rad2deg(atan2(dy, dx) + pi)/361 * 36) + 1; %[-pi, pi]->[0, 2pi]->[0, 360]... -> [1, 36]
            G_magnitude{octave_ind}(:, :, level) = sqrt (dx.^2 + dy.^2);
        elseif ~G_orientation{octave_ind}(:, :, level)
            gaussian_image = G{octave_ind}(:, :, level);
            gaussian_image = padarray(gaussian_image, [sample_radius sample_radius]);
            dx = conv2(gaussian_image, filter_x, 'same'); 
            dy = conv2(gaussian_image, filter_y, 'same');
            G_orientation{octave_ind}(:, :, level) = floor(rad2deg(atan2(dy, dx) + pi)/361 * 36) + 1; %[-pi, pi]->[0, 2pi]->[0, 360]... -> [1, 36]
            G_magnitude{octave_ind}(:, :, level) = sqrt (dx.^2 + dy.^2);
        end

        %%
        orientation_map = G_orientation{octave_ind}(:, :, level);
        magnitude_map = G_magnitude{octave_ind}(:, :, level);
        keypoint_ind = sample_radius + keypoints(ind, 1:2)/2.^(octave_ind-1);
        dominant_orientations = getDominantOrientations(keypoint_ind, sigma, orientation_map, magnitude_map);
        global_descriptors = [global_descriptors; buildDescriptors(keypoint_ind, dominant_orientations, ...
            orientation_map, magnitude_map, sample_radius, cell_size, exp_scale, gauss_mag, desc_array_size)];
        
        
    end
end

function descriptors=buildDescriptors(keypoint_ind, dominant_orientations, ...
    orientation_map, magnitude_map, sample_radius, cell_size, exp_scale, ...
    gauss_mag, desc_array_size)

    descriptors = [];
    start_index_j = keypoint_ind(1) - sample_radius;
    start_index_i = keypoint_ind(2) - sample_radius;
    for i=1:length(dominant_orientations)
        interpdescriptor = zeros(desc_array_size, desc_array_size, 8);
        cos_theta = cosd(dominant_orientations(i));
        sin_theta = sind(dominant_orientations(i));
        rot_mat = [cos_theta -sin_theta;
                   sin_theta cos_theta];
        for y=0:desc_array_size-1
            for x=0:desc_array_size-1
                offset_y = y*cell_size + start_index_j;
                offset_x = x*cell_size + start_index_i;
                for cell_j=offset_y:offset_y+cell_size 
                    for cell_i=offset_x:offset_x+cell_size
                        cell_coord = [cell_j; cell_i] -  keypoint_ind.';
                        cell_coord = rot_mat.' * cell_coord;
                        cell_coord = cell_coord + keypoint_ind.';
                        cell_coord = (cell_coord - [start_index_j; start_index_i])/desc_array_size;

                        if cell_coord(1) <  0 || cell_coord(1) > desc_array_size
                            break;
                        elseif cell_coord(2) < 0 || cell_coord(2) > desc_array_size
                            continue;
                        end
                        dist_j = cell_j - keypoint_ind(1);
                        dist_i = cell_i - keypoint_ind(2);
                        gauss_weight = gauss_mag * exp((dist_j.^2 + dist_i.^2)* exp_scale);
                        pixel_orientation = (interpVal2D(orientation_map, cell_j, cell_i) -1)/35 * 360; %(orientation_map(cell_j, cell_i) -1)/35 * 360;
                        if dist_j == 0 && dist_i == 0
                            shifted_orientation = pixel_orientation;
                        else
                            shifted_orientation = pixel_orientation- dominant_orientations(i);
                        end
                        orientation_index = mod(shifted_orientation, 361)/360 * 8; %[0 8]
                        bin_value =  (interpVal2D(magnitude_map, cell_j, cell_i) * gauss_weight); %(magnitude_map(cell_j, cell_i) * gauss_weight);
                        interpdescriptor = interpolateDescriptor(interpdescriptor, cell_coord, bin_value, orientation_index);
                    end
                end
            end
        end
        %interpdescriptor = reshape(interpdescriptor.', 1, []);
        interpdescriptor = flattenDescriptor(interpdescriptor);
        interpdescriptor = interpdescriptor/norm(interpdescriptor);
        interpdescriptor(interpdescriptor>0.2) = 0.2;
        interpdescriptor = interpdescriptor/norm(interpdescriptor);
        descriptors = [descriptors;interpdescriptor];
    end
    descriptors = unique(descriptors,'rows');
end

function descriptor_flat=flattenDescriptor(descriptor)
    descriptor_flat =[];
    height = size(descriptor, 1);
    width = size(descriptor, 2);
    num_bins = size(descriptor, 3);
    for j=1:height
        for i=1:width
            for bin_index=1:num_bins
                descriptor_flat = [descriptor_flat descriptor(j, i, bin_index)];
            end
        end
    end
end

function descriptor_interp=interpolateDescriptor(descriptor, cell_coord, value, orientation_index)
    height = size(descriptor, 1);
    width = size(descriptor, 2);
    num_bins = size(descriptor, 3);
    descriptor_interp = descriptor;
    for j=1:height
        y_index = abs(cell_coord(1) - (j-0.5));
        if y_index > 1
            continue
        end
        for i=1:width
            x_index = abs(cell_coord(2) - (i-0.5));
            if x_index > 1
                continue
            end
            for o=1:num_bins
                o_index = abs(orientation_index - (o-0.5));
                if o_index > 1
                    continue
                end
                cost = sqrt(y_index.^2 + x_index.^2 + o_index.^2);
                weight = 1 - cost/sqrt(3);
                descriptor_interp(j, i, o) = descriptor_interp(j, i, o) + weight * value;
            end
        end
    end
end

%%get dominant orientations
function dominant_orientations = getDominantOrientations(keypoint_ind, sigma, orientation_candidates, magnitude_map)
    j = keypoint_ind(1);
    i = keypoint_ind(2);
    radius = floor(1.5 * sigma);
    bin = zeros(1, 36);
    exp_scale =  - 1/ (2*sigma.^ 2);
    image_height = size(orientation_candidates, 1);
    image_width = size(orientation_candidates, 2);
    for offset_j=-radius:radius
        y = j + offset_j;
        if y < 1 || y > image_height
            continue;
        end
        for offset_i=-radius:radius
            x = i + offset_i;
            if x < 1 || x > image_width
                continue;
            end
            gauss_weight = exp((offset_j.^2 + offset_i.^2)* exp_scale);
            bin_index = round(interpVal2D(orientation_candidates, y, x)); %orientation_candidates(y, x);
            bin(bin_index) = bin(bin_index) + interpVal2D(magnitude_map, y, x) * gauss_weight; %+ (magnitude_map(y, x) * gauss_weight);
        end
    end
    highest_peak = max(bin);
    highest_peak_ind = find(bin == highest_peak);
    rest_of_bin = bin;
    rest_of_bin(highest_peak_ind) = -1;
    alt_peaks_indices = find(rest_of_bin > (0.8 * highest_peak));
    orientations = [highest_peak_ind alt_peaks_indices];
    dominant_orientations = [];
    for i = 1:length(orientations)
        left_adjacent_orientation = mod(orientations(i) - 2, 36)+1;
        right_adjacent_orientation = mod(orientations(i) , 36) + 1;
        [~, interp_orientation_amount] = interpolatePeak(bin(left_adjacent_orientation), bin(orientations(i)), bin(right_adjacent_orientation));
        orientation = mod(orientations(i) + interp_orientation_amount-1, 35) + 1;
        orientation = ((orientation-1)/35) * 360;
        dominant_orientations = [dominant_orientations orientation];
    end
end

function [peak, peak_ind] = interpolatePeak(left_adjacent, center, right_adjacent)
    peak_ind = 0.5*(left_adjacent - right_adjacent) / (left_adjacent - 2*center + right_adjacent);
    peak = center - 0.25*(left_adjacent - right_adjacent)*peak_ind;
end

function extrema=local_extrema_detection(interval_ind, D_octave_images, octave_ind)
    interval = D_octave_images(:, :, interval_ind:interval_ind+2); %level = interval_ind +1
    image_height = size(interval, 1);
    image_width = size(interval, 2);
    extrema = [];
    for j=2:image_height-1
        for i=2:image_width-1
            pixel_val = interval(j, i, 2);
            is_lowest = false;
            is_highest = false;
            is_equal = false;
            isExtremum = true;
            for level=1:3
                for y=j-1:j+1
                    for x=i-1:i+1
                        if (level == 2) && (y == j) && (x == i)
                            continue
                        end
                        neighbour_pixel_val = interval(y, x, level);
                        if pixel_val < neighbour_pixel_val
                            is_lowest = true;
                        elseif pixel_val > neighbour_pixel_val
                            is_highest = true;
                        elseif pixel_val == neighbour_pixel_val
                            is_equal = true;
                        end
                        if (is_lowest && is_highest) || is_equal
                            isExtremum = false;
                            break
                        end
                    end
                    if ~isExtremum
                        break
                    end
                end
                if ~isExtremum
                    break
                end
            end
            if isExtremum
                pixel_level = interval_ind +1;
                [is_localized, integer_j, integer_i, integer_level, x_hat, hess] = keyPointLocalization(j, i, pixel_level, D_octave_images);
                if is_localized
                    isEliminated = edgeResponseElimination(hess);
                    if ~isEliminated
                        localized_j = integer_j + x_hat(1);
                        localized_i = integer_i + x_hat(2);
                        localized_level = integer_level + x_hat(3);
                        extrema = [extrema; localized_j*2.^(octave_ind-1) localized_i*2.^(octave_ind-1) localized_level octave_ind];
                    end
                end
            end
        end
    end
end

function [is_localized, j, i, level, x_hat, hess]=keyPointLocalization(j, i, level, D_octave_images)
    image_height = size(D_octave_images, 1);
    image_width = size(D_octave_images, 2);
    num_levels = size(D_octave_images, 3);
    is_localized = false;
    pixel_val = D_octave_images(j, i, level);
    grad = [];
    x_hat = 0;
    for step=1:10
        grad = gradient(j, i, level, D_octave_images);
        hess = hessian(j, i, level, D_octave_images);
        if rank(hess) < size(hess, 1) %singular Hessian
            break
        end
        x_hat = -1* (hess\grad);
        if abs(x_hat(1)) <= 0.5 && abs(x_hat(2)) <= 0.5 && abs(x_hat(3)) <= 0.5
            is_localized = true;
            break
        else
            j = j + floor(x_hat(1));
            i = i + floor(x_hat(2));
            level = level + floor(x_hat(3));
            if j < 2 || j > image_height-1 || i < 2 || i > image_width-1 || level < 2 || level > num_levels-1
                break
            end
        end
    end
    if is_localized
        pixel_val = pixel_val + 0.5 * grad.' * x_hat;
        if abs(pixel_val) < 0.03
            is_localized = false;
        end
    end
end

function isEliminated=edgeResponseElimination(hess)
        isEliminated = true;
        r = 10; %threshold
        dxx = hess(1, 1);
        dxy = hess(1, 2);
        dyy = hess(2, 2);
        trace_of_hess = dxx + dyy;
        det_of_hess = dxx*dyy - dxy.^2;
        if (det_of_hess > 0) && (r*trace_of_hess.^2 < det_of_hess * (r+1).^2)
            isEliminated = false;
        end
end

function grad=gradient(j, i, level, D_octave_images)
    dx = (D_octave_images(j, i + 1, level) - D_octave_images(j, i-1, level));
    dy = (D_octave_images(j+1, i, level) - D_octave_images(j-1, i, level));
    dl = (D_octave_images(j, i, level+1) - D_octave_images(j, i, level-1));

    grad = [dx; dy; dl];
end

function hess=hessian(j, i, level, D_octave_images)
    dxx = (D_octave_images(j, i+1, level) - 2*D_octave_images(j, i, level) + D_octave_images(j, i-1, level));
    dyy = (D_octave_images(j+1, i, level) - 2*D_octave_images(j, i, level) +  D_octave_images(j-1, i, level));
    dll = (D_octave_images(j, i, level+1) - 2*D_octave_images(j, i, level) +  D_octave_images(j, i, level-1));
    dxy = (D_octave_images(j+1, i+1, level) - D_octave_images(j-1, i+1, level) - D_octave_images(j+1, i-1, level) + D_octave_images(j-1, i-1, level));
    dxl = (D_octave_images(j, i+1, level+1) - D_octave_images(j, i+1, level-1) - D_octave_images(j, i-1, level+1) + D_octave_images(j, i-1, level-1));
    dyl = (D_octave_images(j+1, i, level+1) - D_octave_images(j-1, i, level+1) - D_octave_images(j+1, i, level-1) + D_octave_images(j-1, i, level-1));

    hess = [dxx dxy dxl;
            dxy dyy dyl;
            dxl dyl dll];
end

function val=interpVal2D(arr, j, i)
    int_j = floor(j);
    int_i = floor(i);
    dec_j = j - int_j;
    dec_i = i - int_i;

    val = (1-dec_j)*(1-dec_i) * arr(int_j, int_i) + ...
            dec_j * (1-dec_i) * arr(int_j+1, int_i) + ...
            (1-dec_j) * dec_i * arr(int_j, int_i +1) + ...
            dec_i * dec_j * arr(int_j + 1, int_i + 1);
end