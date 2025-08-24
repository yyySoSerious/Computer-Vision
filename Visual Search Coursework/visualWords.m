close all;
clear all;

DATASET_FOLDER = '../MSRC_ObjCategImageDatabase_v2';
DESCRIPTOR_FOLDER = '../descriptors';
DESCRIPTOR_SUBFOLDER='sift';
OUT_FOLDER = '../descriptors';
OUT_SUBFOLDER='visual_words';
folder = [OUT_FOLDER,'/',OUT_SUBFOLDER,'/'];
if ~exist(folder, 'dir')
    mkdir(folder);
end

ALLFILES=cell(1,0);
ctr=1;
allfiles=dir (fullfile([DATASET_FOLDER,'/Images/*.bmp']));
num_images = length(allfiles);
num_descriptors=zeros(num_images, 1);
local_descriptors=[];
fprintf("\nLoading descriptors...");
for filenum=1:num_images
    fname=allfiles(filenum).name;
    imgfname_full=([DATASET_FOLDER,'/Images/',fname]);
    featfile=[DESCRIPTOR_FOLDER,'/',DESCRIPTOR_SUBFOLDER,'/',fname(1:end-4),'.mat'];%replace .bmp with .mat
    load(featfile,'F');
    ALLFILES{ctr}=imgfname_full;
    num_descriptors(ctr) = size(F, 1);
    local_descriptors=[local_descriptors ; F]; %size: n  x ncolsof('F')
    ctr=ctr+1;
end

max_iterations = floor(0.1 * num_images);
fprintf("\nComputing Kmeans...");

[visual_words, num_occurences] = computeKMeans(local_descriptors, 60, 1000); 
[visual_words, num_occurences] = removeCommonWords(visual_words, num_occurences, 0.02);

inv_doc_freq = log10(repmat(sum(num_occurences), 1, length(num_occurences)) ./num_occurences);
fprintf('\nSaving descriptors...\n');
for i=1:num_images
    if i ==1
        start_index = 1;
        end_index = num_descriptors(i);
    else
        start_index = sum(num_descriptors(1:i-1))+1;
        end_index = start_index + num_descriptors(i)-1;
    end

    fname=allfiles(i).name;
    fprintf('Saving descriptor %d/%d - %s\n',i,length(allfiles),fname);
    tic;
    imgfname_full=([DATASET_FOLDER,'/Images/',fname]);
    img=double(imread(imgfname_full));
    fout=[folder, fname(1:end-4),'.mat'];%replace .bmp with .mat
    F=createGlobalDescriptor(local_descriptors(start_index:end_index, :), visual_words, inv_doc_freq);
    save(fout,'F');
    toc
end

function global_descriptor=createGlobalDescriptor(local_descriptors, visual_words, inv_doc_freq)
    num_local_descriptors = size(local_descriptors, 1);
    num_words = size(visual_words, 1);
    global_descriptor = zeros(1, num_words);
    for ind=1:num_local_descriptors
            dist = visual_words(: , :) - repmat(local_descriptors(ind,:), size(visual_words, 1), 1);
            dist = vecnorm(dist, 2, 2);
            [~, word_ind] = min(dist);
            global_descriptor(word_ind) = global_descriptor(word_ind) + 1;
    end
    non_empty_bins = global_descriptor(global_descriptor>0);
    global_descriptor = (global_descriptor ./ length(non_empty_bins)) .* inv_doc_freq;
    global_descriptor = global_descriptor ./ norm(global_descriptor);
end

function [visual_words, num_occurences]=removeCommonWords(visual_words, num_occurences, percentage)
    num_words_to_remove = round(percentage * size(visual_words, 1));
    for i=1:num_words_to_remove
        if length(num_occurences) < 1
            break;
        end
        [~, ind] = max(num_occurences);
        if ind > 1 && ind < size(visual_words, 1)
            visual_words= [visual_words(1: ind-1, :);visual_words(ind+1:end, :)];
            num_occurences= [num_occurences(1: ind-1) num_occurences(ind+1:end)];
        elseif ind == 1
            visual_words = visual_words(2:end, :);
            num_occurences = num_occurences(2:end);
        elseif ind == size(visual_words, 1)
            visual_words = visual_words(1:end-1, :);
            num_occurences = num_occurences(1:end-1);
        end
    end
end