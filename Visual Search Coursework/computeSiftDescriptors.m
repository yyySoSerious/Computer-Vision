close all;
clear all;

DATASET_FOLDER = '../MSRC_ObjCategImageDatabase_v2';
OUT_FOLDER = '../descriptors';
OUT_SUBFOLDER='sift';

folder = [OUT_FOLDER,'/',OUT_SUBFOLDER,'/'];
if ~exist(folder, 'dir')
    mkdir(folder);
end
allfiles=dir (fullfile([DATASET_FOLDER,'/Images/*.bmp']));
num_intervals = 2;
num_octaves = -1;
for filenum=1:length(allfiles)
    fname=allfiles(filenum).name;
    fprintf('Processing file %d/%d - %s\n',filenum,length(allfiles),fname);
    tic;
    imgfname_full=([DATASET_FOLDER,'/Images/',fname]);
    img=imread(imgfname_full);
    fout=[folder, fname(1:end-4),'.mat'];%replace .bmp with .mat
    F=siftDescriptors(img, num_intervals, num_octaves);
    save(fout,'F');
    toc
end
