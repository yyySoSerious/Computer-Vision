function [descs_lowerDim, eig_vals]=pca_reduced(original_descriptors, numSignificants)
    num_descriptors  = size(original_descriptors, 2);
    mu = mean(original_descriptors, 2);
    descriptors_sub = original_descriptors - repmat(mu, 1, num_descriptors);
    covarr_mat = (descriptors_sub*descriptors_sub') ./num_descriptors;
    [V, D] = eig(covarr_mat);
    eig_vals = diag(D);
    V_reduced = [];
    if isempty(numSignificants)
        trace = sum(eig_vals);
        significant_trace = 0.97 * trace;
        significant_ind = [];
        total = 0;
        for i=length(eig_vals):-1:1
            significant_ind = [significant_ind i];
            total = total + eig_vals(i);
            if total > significant_trace
                break;
            end
        end
        significant_ind = flip(significant_ind);
        V_reduced = V(:, significant_ind);
        eig_vals = eig_vals(significant_ind);
    else
        V_reduced = V(:, end-numSignificants+1:end);
        eig_vals = eig_vals(end-numSignificants+1:end);
    end

    %% project descriptors onto a lower dimensional space
    descs_lowerDim = (V_reduced' * original_descriptors);
end