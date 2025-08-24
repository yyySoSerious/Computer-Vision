%% EEE3032 - Computer Vision and Pattern Recognition (ee3.cvpr)
%%
%% cvpr_computedescriptors.m
%% Skeleton code provided as part of the coursework assessment
%% This code will iterate through every image in the MSRCv2 dataset
%% and call a function 'extractRandom' to extract a descriptor from the
%% image.  Currently that function returns just a random vector so should
%% be changed as part of the coursework exercise.
%%
%% (c) John Collomosse 2010  (J.Collomosse@surrey.ac.uk)
%% Centre for Vision Speech and Signal Processing (CVSSP)
%% University of Surrey, United Kingdom

close all;
clear all;

%% Edit the following line to the folder you unzipped the MSRCv2 dataset to
DATASET_FOLDER = '../MSRC_ObjCategImageDatabase_v2';

%% Create a folder to hold the results...
OUT_FOLDER = '../descriptors';
%% and within that folder, create another folder to hold these descriptors
%% the idea is all your descriptors are in individual folders - within
%% the folder specified as 'OUT_FOLDER'.
descriptor_types = ["globalRGBhisto", "spatialGrid", "eoh", "eohWithColor"];
OUT_SUBFOLDER = descriptor_types(1);
%OUT_SUBFOLDER = descriptor_types(2);
%OUT_SUBFOLDER = descriptor_types(3);
%OUT_SUBFOLDER = descriptor_types(4);

OUT_SUBFOLDER = char(OUT_SUBFOLDER);
folder = [OUT_FOLDER,'/',OUT_SUBFOLDER,'/'];
if ~exist(folder, 'dir')
    mkdir(folder)
end
allfiles=dir (fullfile([DATASET_FOLDER,'/Images/*.bmp']));
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
    save(fout,'F');
    toc
end
