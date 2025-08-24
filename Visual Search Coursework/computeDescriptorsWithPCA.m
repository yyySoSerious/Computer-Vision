close all;
clear all;

DATASET_FOLDER = '../MSRC_ObjCategImageDatabase_v2';
OUT_FOLDER = '../descriptors';
descriptor_types = ["globalRGBhisto", "spatialGrid", "eoh", "eohWithColor"];
%OUT_SUBFOLDER = descriptor_types(1);
%OUT_SUBFOLDER = descriptor_types(2);
%OUT_SUBFOLDER = descriptor_types(3);
OUT_SUBFOLDER = descriptor_types(4);

OUT_SUBFOLDER = char(OUT_SUBFOLDER);
folder = [OUT_FOLDER,'/',OUT_SUBFOLDER,'/'];
if ~exist(folder, 'dir')
    mkdir(folder)
end
allfiles=dir (fullfile([DATASET_FOLDER,'/Images/*.bmp']));
descriptors = [];
num_classes = 20;
library_classes = cell(num_classes, 1);
for filenum=1:length(allfiles)
    fname=allfiles(filenum).name;
    fprintf('Processing file %d/%d - %s\n',filenum,length(allfiles),fname);
    tic;

    imgfname_full=([DATASET_FOLDER,'/Images/',fname]);
    img=double(imread(imgfname_full));
    fout=[folder, fname(1:end-4),'.mat'];%replace .bmp with .mat
    switch OUT_SUBFOLDER
        case descriptor_types(1)
            F = globalColorHistogram(img, 5); %extractRandom(img);
        case descriptor_types(2)
            F = spatialGrid(img, 5);
        case descriptor_types(3)
            F = eoh(img, 5, 18, 0.06);
        case descriptor_types(4)
            F = eohWithColor(img, 5, 17, 0.06);
        otherwise
            F = extractRandom(img);
    end
    descriptors = [descriptors F'];
    toc
end

%% project descriptors to 3-dimensional space and save relevant eigen values
%[descriptors, eig_vals] = pca_reduced(descriptors, 3); 

%% project descriptors onto lower dimensional space, which is determined by the trace of D
[descriptors, eig_vals] = pca_reduced(descriptors, []); 

%% save eigen values
eig_valsout = [folder, 'eig_vals.mat'];
save(eig_valsout, 'eig_vals');

%% save lower dimensional descriptors onto file and group the descriptors according to their class
for filenum=1:length(allfiles)
        fname=allfiles(filenum).name;
        fprintf('Processing file %d/%d - %s\n',filenum,length(allfiles),fname);
        tic;
        fout=[folder, fname(1:end-4),'.mat'];%replace .bmp with .mat
        F = descriptors(:, filenum)';
        save(fout,'F');
        
        someClass = getfield(split(fname, '_'), {1});
        className = str2double(someClass{1});
        library_classes{className, 1} = [library_classes{className, 1} F'];
        toc
end

%% use the grouped-up descriptors to generate eigen models for the class and save the models 
eigModels = buildEigenModels(library_classes); %an eigModel = [V, D, mu]' 
eigModelsout = [folder, 'eigModels.mat'];
save(eigModelsout, 'eigModels');
