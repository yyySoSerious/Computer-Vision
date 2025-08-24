function eigModels=buildEigenModels(library_classes)
    num_classes = size(library_classes, 1);
    eigModels = cell(20, 1);
    for class_ind=1:num_classes
        classDescriptors = library_classes{class_ind, 1};
        num_descriptors  = size(classDescriptors, 2);
        mu = mean(classDescriptors, 2);
        classDescriptors_sub = classDescriptors - repmat(mu, 1, num_descriptors);
        covarr_mat = (classDescriptors_sub*classDescriptors_sub') ./num_descriptors;
        [V, D] = eig(covarr_mat);
        eigModels{class_ind, 1} = {V;D;mu};
    end
end