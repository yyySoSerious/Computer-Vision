function [centroids, numElementsInCentroid] = computeKMeans(descriptors, k, max_iterations)
    num_descriptors = size(descriptors, 1);
    descriptor_dim = size(descriptors, 2);

    randInd = randperm(num_descriptors);
    k = min([length(randInd), k]);
    chosen_descriptors = descriptors;
    chosen_descriptors(randInd, :) = descriptors;
    centroids = chosen_descriptors(1:k, :);
    counter = 0;
    while true
        clusters = cell(1, k);
        
        %move descriptors to their nearest cluster
        tic
        for ind=1:num_descriptors
            dist = centroids - repmat(descriptors(ind,:), size(centroids, 1), 1);
            dist = vecnorm(dist, 2, 2);
            [~, closest_cluster_ind] = min(vecnorm(dist, 2, 2));
            clusters{closest_cluster_ind} = [clusters{closest_cluster_ind}; descriptors(ind, :)];
        end
        toc;

        counter = counter + 1;
        old_centroids = centroids;
        %adjust centroid for each cluster
        for ind=1:k
            centroids(ind, :) = mean(clusters{ind}, 1);
        end
        diff = sum(old_centroids - centroids, "all");

        if ~diff || counter > max_iterations
            %count number of elements in cluster before exiting
            numElementsInCentroid = zeros(1, k);
            for ind=1:k
                numElementsInCentroid(ind) = size(clusters{ind}, 1);
            end
            fprintf("\nnum iterations: %d", counter);
            break;
        end
    end
end