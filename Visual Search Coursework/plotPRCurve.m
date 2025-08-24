function plotPRCurve(FeatDesc, filenames, distMetric, eigModels, eig_vals)
    numDescriptors = size(FeatDesc, 1);
    randInd = randperm(numDescriptors);
    shuffledFeatDescriptors = zeros(size(FeatDesc));
    shuffedFilenames = cell(size(filenames));
    shuffledFeatDescriptors(randInd, :) = FeatDesc;
    shuffedFilenames(randInd) = filenames;
    queryRatio = 0.1;
    nQueryInd = max(1e-8, floor(queryRatio * numDescriptors));
    nQueries = floor(numDescriptors / nQueryInd);
    library_sz = numDescriptors - nQueryInd;
    avg_precision = zeros(1, library_sz+1);
    avg_recall = zeros(1, library_sz+1);
    avg_map = 0;
    visualParams = {};%Parameters for visualize function
    for queryblock=0:nQueries-1
        startQueryIndex = nQueryInd*queryblock + 1;
        endQueryIndex = nQueryInd*queryblock + nQueryInd;

        queries = shuffledFeatDescriptors(startQueryIndex:endQueryIndex, :);
        queries_filenames = shuffedFilenames(startQueryIndex:endQueryIndex);
        
        if queryblock == 0
            library = shuffledFeatDescriptors(endQueryIndex +1:numDescriptors, :);
            library_filenames = shuffedFilenames(endQueryIndex +1:numDescriptors);

        elseif endQueryIndex < numDescriptors
            library = [shuffledFeatDescriptors(1:startQueryIndex-1, :); 
                shuffledFeatDescriptors(endQueryIndex + 1:numDescriptors, :)];
            library_filenames = {shuffedFilenames{1:startQueryIndex-1} shuffedFilenames{endQueryIndex + 1:numDescriptors}};
        else
            library = shuffledFeatDescriptors(1:startQueryIndex-1, :);
            library_filenames = shuffedFilenames(1:startQueryIndex-1);
        end

        libraryClasses = getLibraryClasses(library_filenames);
        queries_sz = size(queries, 1);

        avg_block_precision = zeros(1, library_sz+1);
        avg_block_recall = zeros(1, library_sz+1);
        block_map = 0;

        fprintf('\n\n------------------------');
        fprintf('\n------------------------');
        fprintf('\nBlOCK');
        fprintf('\n------------------------');
        fprintf('\n------------------------');
        for i=1:queries_sz
            query = queries(i, :);
            query_filename = queries_filenames{i};
            query_class = getfield(split(getfield(split(query_filename, '/'), {4}), '_'), {1});
            query_class = str2double(query_class{1});

            relevant = libraryClasses == query_class;
            relevant_size = length(find(relevant));
            dst = getSimilarity(query, query_class, library, libraryClasses, distMetric, eigModels, eig_vals);
            dst_classes = dst(:, 2);
            query_result_ind = find(dst_classes == 0);
            if query_result_ind == 1
                results = [query_class libraryClasses(dst_classes(2:end))];
            elseif query_result_ind == length(dst_classes)
                results = [libraryClasses(dst_classes(1:end-1)) query_class];
            else
                results = [libraryClasses(dst_classes(1:query_result_ind-1)) query_class libraryClasses(dst_classes(query_result_ind+1:end))];
            end
            %results = libraryClasses(dst(:, 2) ~= 0);

            fprintf('\n\n-----');
            fprintf('\nQuery');
            fprintf('\n-----');
            fprintf('\nquery class: %d', query_class);
            fprintf('\nnumber of relevant data in library: %d', relevant_size);
            fprintf('\nresults of the classes are:    ');
            fprintf('%-6d ', results);
            results = results == query_class;
            fprintf('\nrelevant results are:          ');
            fprintf('%-6d ', results);
            num_correct = [];
            for j=1:length(results)
                num_correct = [num_correct sum(results(1:j))];
            end
            fprintf('\n# of top N that are correct:   ');
            fprintf('%-7d', num_correct);
            top_n = 1:length(results);
            query_precision = num_correct ./ top_n;
            query_recall = num_correct ./ relevant_size;
            query_recall(query_recall > 1) = 1;
            fprintf('\nquery precision:               ');
            fprintf('%-6.3f ', double(query_precision));
            fprintf('\nquery recall:                  ');
            fprintf('%-6.3f ', double(query_recall));
            avg_block_precision = avg_block_precision + query_precision;
            avg_block_recall = avg_block_recall + query_recall;
            fprintf('\ntotal block precision:         ');
            fprintf('%-6.3f ', double(avg_block_precision));
            fprintf('\ntotal block recall:            ');
            fprintf('%-6.3f ', double(avg_block_recall));
            ap_num = query_precision .* results; %numerator for ap
            query_ap = sum(ap_num) / relevant_size;
            block_map = block_map + query_ap;

            %save details for visualizing the results of first query of first block 
            if queryblock == 0 && i == 1
                visualParams = {query, query_filename, query_class, library, library_filenames, libraryClasses, eigModels};
                %visualizeResults(query, query_filename, query_class , library, library_filenames, libraryClasses);
            end
        end
        avg_block_precision = avg_block_precision ./ queries_sz;
        avg_block_recall = avg_block_recall ./ queries_sz;
        fprintf('\n\n-------------');
        fprintf('\nBlock Results');
        fprintf('\n-------------');
        fprintf('\naverage block precision:       ');
        fprintf('%-6.3f ', double(avg_block_precision));
        fprintf('\naverage block recall:          ');
        fprintf('%-6.3f ', double(avg_block_recall));
        block_map = block_map ./queries_sz;
        fprintf('\nblock map:                     ');
        fprintf('%-6.3f', double(block_map));
        avg_precision = avg_precision + avg_block_precision;
        avg_recall = avg_recall + avg_block_recall;
        avg_map = avg_map + block_map;

    end
    avg_precision = avg_precision ./ nQueries;
    avg_recall = avg_recall ./ nQueries;
    avg_map = avg_map ./ nQueries;

    fprintf('\n\n------------------------');
    fprintf('\n------------------------');
    fprintf('\nFINAL RESULTS');
    fprintf('\n------------------------');
    fprintf('\n------------------------');
    fprintf('\naverage precision:             ');
    fprintf('%-6.3f ', double(avg_precision));
    fprintf('\naverage recall:                ');
    fprintf('%-6.3f ', double(avg_recall));
    fprintf('\naverage map:                   ');
    fprintf('%-6.3f ', double(avg_map));
    plot(avg_recall, avg_precision);
    xlabel('Recall');
    ylabel('Precision');
    axis([0 1 0 1])

    %display results of the first query of the first block
    visualizeResults(visualParams{1}, visualParams{2}, visualParams{3}, ...
        visualParams{4}, visualParams{5}, visualParams{6}, distMetric, visualParams{7}, eig_vals);

end

function img=readEditImg(filename)
   img=imread(filename);
   img=img(1:2:end,1:2:end,:); % make image a quarter size
   img=img(1:81,:,:); % crop image to uniform size vertically (some MSVC images are different heights)
end

function dst=getSimilarity(query, query_class, library, libraryClasses, distMetric, eigModels, eig_vals)
    library_sz = size(library, 1);
    eigModel = [];
    if not(isempty(eigModels))
        eigModel = eigModels{query_class, 1}; %candidate's eigen model
    end
    dst = [];
    dst = [dst; [cvpr_compare(query, query, distMetric, eigModel, eig_vals) 0]];
    for i=1:library_sz
        candidate = library(i,:);
        candidateClass = libraryClasses(i);
        eigModel = [];
        if not(isempty(eigModels))
            eigModel = eigModels{candidateClass, 1}; %candidate's eigen model
        end
        thedst=cvpr_compare(query, candidate, distMetric, eigModel, eig_vals);
        dst=[dst ; [thedst i]];
    end
    dst=sortrows(dst,1);  % sort the results
end

function visualizeResults(query, queries_filename, query_class, library, library_filenames, libraryClasses, distMetric, eigModels, eig_vals)
    dst = getSimilarity(query, query_class, library,  libraryClasses, distMetric, eigModels, eig_vals);
    fprintf('\n\nVisualizing a query belonging to class %d\n', query_class);
    queryImg = readEditImg(queries_filename);
    SHOW=15; % Show top 15 results
    dst=dst(1:SHOW,:);
    topnDisplay=[];
    for i=1:size(dst,1)
        class_id = dst(i,2);
        if class_id == 0
            img = queryImg;
        else
            img = readEditImg(library_filenames{dst(i,2)});
        end
        topnDisplay=[topnDisplay img];
    end
    figure;
    ax1 = axes('Position',[0.3 0.5 0.4 0.4]);
    ax2 = axes('Position',[0.01 0.01 0.99 0.5]);

    axes(ax1);
    imshow(queryImg);
    axes(ax2);
    imshow(topnDisplay);%montage({topnDisplay},"Size",[3 5]);
    axis off;
end

function libraryClasses=getLibraryClasses(library_filenames)
    library_sz = size(library_filenames, 2);
    libraryClasses = zeros(1, library_sz);
    for i=1:library_sz
        someClass = getfield(split(getfield(split(library_filenames{i}, '/'), {4}), '_'), {1});
        libraryClasses(i) = str2double(someClass{1});
    end
end
