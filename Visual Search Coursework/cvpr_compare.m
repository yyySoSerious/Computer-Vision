function dst=cvpr_compare(F1, F2, distMetric, eigModel, eig_vals)

% This function should compare F1 to F2 - i.e. compute the distance
% between the two descriptors

% For now it just returns a random number

    %distMetric: EUCLIDEAN, MAHALANOBIS

    if strcmp(distMetric, 'EUCLIDEAN')
        dst = euclidean_distance(F1, F2);
    elseif strcmp(distMetric, 'L1NORM')
        dst = l1_distance(F1, F2);
    elseif strcmp(distMetric, 'BHATTACHARYYA')
        dst = bhattacharyya_distance(F1, F2);
    elseif strcmp(distMetric, 'JEFFRIES-MATUSITA')
        dst = jefMat_distance(F1, F2);
    elseif strcmp(distMetric, 'COSINE')
        dst = cosine_distance(F1, F2);
    elseif strcmp(distMetric, 'MAHALANOBIS')
        dst = mahalanobis_distance(F1, F2, eig_vals);
    elseif strcmp(distMetric, 'costMAHALANOBIS')
        dst = mahalanobis_distanceCost(F1, F2, eigModel);
    end

end


function dst=euclidean_distance(F1, F2)
    x = F1 - F2;
    x = x.^2;
    x = sum(x);
    dst = sqrt(x);
end

function dst=l1_distance(F1, F2)
    x = abs(F1 - F2);
    dst = sum(x);
end

function dst=bhattacharyya_distance(F1, F2)
    unitF1 = F1 / sqrt(F1 * F1');
    unitF2 = F2 / sqrt(F2 * F2');
    x = sqrt(unitF1 .* unitF2);
    x = sum(x);
    dst = -log(x);
end

%Jeffries-Matusita Distance
function dst=jefMat_distance(F1, F2)
    x = 2 * (1 - exp(-bhattacharyya_distance(F1, F2)));
    dst = sqrt(x);
end

function dst=cosine_distance(F1, F2)
    normF1 = sqrt(F1 * F1');
    normF2 = sqrt(F2 * F2');
    dot_x = F1 * F2';
    dst = -dot_x / (normF1 * normF2);
end

function dst=mahalanobis_distance(F1, F2, eig_vals)
    x = F1 - F2;
    x = x.^2;
    x = x ./ flip(eig_vals');
    x = sum(x);
    dst = sqrt(x);
end
function dst=mahalanobis_distanceCost(F1, F2, eigModel)
    F1 = F1';
    F2 = F2';
    V = eigModel{1};
    D = eigModel{2};
    mu = eigModel{3};
    F1_sub = F1 - mu;
    cost = sqrt(F1_sub' * V * diag(1 ./ diag(D)) * V' * F1_sub);
    dst = euclidean_distance(F1, F2)*cost;
end

