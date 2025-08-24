%% EEE3032 - Computer Vision and Pattern Recognition (ee3.cvpr)
%%
%% cvpr_visualsearch.m
%% Skeleton code provided as part of the coursework assessment
%%
%% This code will load in all descriptors pre-computed (by the
%% function cvpr_computedescriptors) from the images in the MSRCv2 dataset.
%%
%% It will pick a descriptor at random and compare all other descriptors to
%% it - by calling cvpr_compare.  In doing so it will rank the images by
%% similarity to the randomly picked descriptor.  Note that initially the
%% function cvpr_compare returns a random number - you need to code it
%% so that it returns the Euclidean distance or some other distance metric
%% between the two descriptors it is passed.
%%
%% (c) John Collomosse 2010  (J.Collomosse@surrey.ac.uk)
%% Centre for Vision Speech and Signal Processing (CVSSP)
%% University of Surrey, United Kingdom

close all;
clear all;

%% Edit the following line to the folder you unzipped the MSRCv2 dataset to
DATASET_FOLDER = '../MSRC_ObjCategImageDatabase_v2';

%% Folder that holds the results...
DESCRIPTOR_FOLDER = '../descriptors';
%% and within that folder, another folder to hold the descriptors
%% we are interested in working with
DESCRIPTOR_SUBFOLDER='globalRGBhisto';
%DESCRIPTOR_SUBFOLDER='spatialGrid';
%DESCRIPTOR_SUBFOLDER='eoh';
%DESCRIPTOR_SUBFOLDER='eohWithColor';

%DESCRIPTOR_SUBFOLDER='visual_words';


%% 1) Load all the descriptors into "ALLFEAT"
%% each row of ALLFEAT is a descriptor (is an image)

ALLFEAT=[];
ALLFILES=cell(1,0);
ctr=1;
allfiles=dir (fullfile([DATASET_FOLDER,'/Images/*.bmp']));
for filenum=1:length(allfiles)
    fname=allfiles(filenum).name;
    imgfname_full=([DATASET_FOLDER,'/Images/',fname]);
    featfile=[DESCRIPTOR_FOLDER,'/',DESCRIPTOR_SUBFOLDER,'/',fname(1:end-4),'.mat'];%replace .bmp with .mat
    load(featfile,'F');
    ALLFILES{ctr}=imgfname_full;
    ALLFEAT=[ALLFEAT ; F]; %size: n  x ncolsof('F')
    ctr=ctr+1;
end

%% plot PR Curve and visualize the first query of cross-validation test,
%% using different distance measures
plotPRCurve(ALLFEAT, ALLFILES, 'EUCLIDEAN', [], []);
%plotPRCurve(ALLFEAT, ALLFILES, 'L1NORM', [], []);
%plotPRCurve(ALLFEAT, ALLFILES, 'BHATTACHARYYA', [], []);
%plotPRCurve(ALLFEAT, ALLFILES, 'JEFFRIES-MATUSITA', [], []);
%plotPRCurve(ALLFEAT, ALLFILES, 'COSINE', [], []);

%% Mahalanobis distance measure
%eig_vals_file = [DESCRIPTOR_FOLDER,'/',DESCRIPTOR_SUBFOLDER,'/','eig_vals.mat'];
%load(eig_vals_file, 'eig_vals');
%plotPRCurve(ALLFEAT, ALLFILES, 'MAHALANOBIS', [], eig_vals);

%% Experimenting with using Mahalanobis Distance as a cost metric
%% computeDescriptorsWithPCA should be run first if the eigen models for
%% each class has not been saved.
%eigModels_file = [DESCRIPTOR_FOLDER,'/',DESCRIPTOR_SUBFOLDER,'/','eigModels.mat'];
%load(eigModels_file, 'eigModels');
%plotPRCurve(ALLFEAT, ALLFILES, 'costMAHALANOBIS', eigModels, []); 
