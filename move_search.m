function [match_success, i_end, j_end, s_end, S_NaN, search_steps, search_path] = move_search(im, w, i_start, j_start, search_radius, varargin)
%MOVE_SEARCH Searches a matching point of template, W, in image, R,
%from a approximate matching point that stores in varialbe, LAYER_BEST.
%Pay attention that MOVE_SEARCH can't always find the best matching point.
%Matching will carry out in a searching frame that will move step by step
%closing to target point. If the best matching point is in the centre of
%the searching box, searching end.
%
%   SEARCH_PATH records the starting point and (SEARCH_STEPS) searched
%   points.
% 
%   POLARITY stands for similarity metric polarity of matching. POLARITY
%   can choose from below options, the default option is enclosed in
%   brace({}):
%       {'minimum'}     The smaller the similarity metric, the better the
%                       matching.
%       'maximum'       The bigger the similarity metric, the better the
%                       matching.
% 
% Example: [i_end, j_end, s_end, S_NaN, search_steps, search_path] = move_search(im, w, i_start, j_start, search_radius, varargin) 

% Start timing.
tic;

% Check right number of inputs.
error(nargchk(5, 6, nargin));

% Set default values.
polarity = 'minimum';

% Get inputs.
if nargin > 5
    polarity = varargin{1};
end

[im_m, im_n] = size(im);
match_success = 1;

% Count the searching steps.
search_steps = 0;

% Sign of centreing.
centre = 1;

% Make a similarity matrix filled with NaN.
S_NaN = NaN(size(im));

% Pad image, IM, for convolution operation.
im_padded = pad_image(im, w);

% Keep the search-starting point in search_path.
search_path = {[i_start, j_start]};
    
% When centred, stop searching.
while(centre)
    % Move_step add one.
    search_steps = search_steps + 1;
    
    % Make the searching box.
    row_start = i_start-search_radius;
    row_end = i_start+search_radius;
    column_start = j_start-search_radius;
    column_end = j_start+search_radius;
    
    if row_start < 1 || row_end > im_m || column_start < 1 || column_end > im_n
        warning('MOVE_SEARCH:BeyondBoundary', 'This move_search failed due to the matching target is too close too the boundary')
        match_success = 0;
        i_end = NaN;
        j_end = NaN;
        s_end = NaN;
        
        figure, imshow(S_NaN, [ ]);
        
        text(search_path{1}(1,2), search_path{1}(1,1),'o','color','r');
        % Add one pixel to label point, to keep a small distance between target 
        % point and label text for better display.
        text(search_path{1}(1,2)+0.2, search_path{1}(1,1)+0.2,['Start(0):(',num2str(search_path{1}(1,1)),',',num2str(search_path{1}(1,2)),')'],'color','b');
        for tt = 2:search_steps-1
            text(search_path{tt}(1,2), search_path{tt}(1,1),'o','color','r');
            text(search_path{tt}(1,2)+0.2, search_path{tt}(1,1)+0.2,['(',num2str(tt),')'],'color','b');
        %     text(search_path{tt}(1,2)+1, search_path{tt}(1,1)+1,['(',num2str(tt),'):(',num2str(search_path{tt}(1,2)),',',num2str(search_path{tt}(1,1)),')'],'color','b');
        end
        text(search_path{end}(1,2), search_path{end}(1,1),'o','color','r');
        %??tt??text(search_path{end}(1,2)+0.2, search_path{end}(1,1)+0.2,['End(',num2str(tt+1),'):(',num2str(search_path{end}(1,1)),',',num2str(search_path{end}(1,2)),')'],'color','b');
        
        return
    end
    
    %Fill similarity values.
    for ii2 = row_start : row_end
        for jj2 = column_start : column_end
            
            % Debugging code:
%             if isnan(ii2)
%                 ii2 = NaN;
%             end
            % Debugging code:
            
            if isnan(S_NaN(ii2, jj2))
                S_NaN(ii2, jj2) = point_similarity(im, im_padded, w, ii2, jj2);
            else
                continue
            end
        end
    end
    
    % Get the sub-window image of searching frame from similarity matrix.
    search_box = S_NaN(row_start:row_end, column_start:column_end);

    % Find similarity peak location.
    switch polarity
        case 'minimum'
            s_end = min(search_box(:));
            [i_end, j_end] = find(S_NaN == s_end);
            if max(size(i_end)) ~= 1
                error('The peak-similarity is not unique in the searching box.(move_search)');
            end
        case 'maximum'
            s_end = max(search_box(:));
            [i_end, j_end] = find(S_NaN == s_end);
            if max(size(i_end)) ~= 1
                error('The peak-similarity is not unique in the searching box.(move_search)');
            end            
        otherwise
            error('Nnknown option for POLARITY.(move_search)');
    end
    
    % Keep searched points.
    search_path{search_steps+1} = [i_end, j_end];
    
    % Check whether centreing.   
    if i_end==i_start && j_end==j_start
        centre = 0;
    else
        i_start = i_end;
        j_start = j_end;
    end
    
end

% End timing.
disp('Function: move_search');
toc;

% Distance check.
if max([abs(i_end - i_start), abs(j_end - j_start)]) > min(size(w))
    warning('MOVE_SEARCH:Distance_check','Move_search distance is larger than width of template.There is a risk of matching another target.')
end

% Output.

disp(['Start:(', num2str(search_path{1}(1,1)), ', ', num2str(search_path{1}(1,2)), ') ']);
disp(['End:(', num2str(i_end), ', ', num2str(j_end), ') ']);
disp(['matched-similarity: ',num2str(s_end)]);
disp(['searched-steps: ',num2str(search_steps)]);


figure, imshow(S_NaN, [ ]);

text(search_path{1}(1,2), search_path{1}(1,1),'o','color','r');
% Add one pixel to label point, to keep a small distance between target 
% point and label text for better display.
text(search_path{1}(1,2)+0.2, search_path{1}(1,1)+0.2,['Start(0):(',num2str(search_path{1}(1,1)),',',num2str(search_path{1}(1,2)),')'],'color','b');
for tt = 2:search_steps-1
    text(search_path{tt}(1,2), search_path{tt}(1,1),'o','color','r');
    text(search_path{tt}(1,2)+0.2, search_path{tt}(1,1)+0.2,['(',num2str(tt),')'],'color','b');
%     text(search_path{tt}(1,2)+1, search_path{tt}(1,1)+1,['(',num2str(tt),'):(',num2str(search_path{tt}(1,2)),',',num2str(search_path{tt}(1,1)),')'],'color','b');
end
text(search_path{end}(1,2), search_path{end}(1,1),'o','color','r');
%??tt??text(search_path{end}(1,2)+0.2, search_path{end}(1,1)+0.2,['End(',num2str(tt+1),'):(',num2str(search_path{end}(1,1)),',',num2str(search_path{end}(1,2)),')'],'color','b');



% % Display.
% S_T = S_NaN';
% figure, surf(S_T(1:round(end/50):end, 1:round(end/50):end));

% Convert cell, SEARCH_PATH, to matrix.
search_path = cell2mat(search_path);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time VS. Search_radius %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% [i_end1, j_end1, s_end1, S_NaN1, search_steps1, search_path1] = move_search(f_squr, w2, 104, 150, 1);
% Function: move_search
% Elapsed time is 0.000672 seconds.
% (104, 155) 
% matched-similarity: 240440
% searched-steps: 6

% [i_end1, j_end1, s_end1, S_NaN1, search_steps1, search_path1] = move_search(f_squr, w2, 104, 150, 2);
% Function: move_search
% Elapsed time is 0.000873 seconds.
% (109, 155) 
% matched-similarity: 58813
% searched-steps: 6

% [i_end1, j_end1, s_end1, S_NaN1, search_steps1, search_path1] = move_search(f_squr, w2, 104, 150, 3);
% Function: move_search
% Elapsed time is 0.000823 seconds.
% (109, 155) 
% matched-similarity: 58813
% searched-steps: 3

% [i_end1, j_end1, s_end1, S_NaN1, search_steps1, search_path1] = move_search(f_squr, w2, 104, 150, 4);
% Function: move_search
% Elapsed time is 0.000516 seconds.
% (109, 155) 
% matched-similarity: 58813
% searched-steps: 3

% [i_end1, j_end1, s_end1, S_NaN1, search_steps1, search_path1] = move_search(f_squr, w2, 104, 150, 5);
% Function: move_search
% Elapsed time is 0.000576 seconds.
% (109, 155) 
% matched-similarity: 58813
% searched-steps: 2

% [i_end1, j_end1, s_end1, S_NaN1, search_steps1, search_path1] = move_search(f_squr, w2, 104, 150, 6);
% Function: move_search
% Elapsed time is 0.001174 seconds.
% (109, 155) 
% matched-similarity: 58813
% searched-steps: 2

% [i_end1, j_end1, s_end1, S_NaN1, search_steps1, search_path1] = move_search(f_squr, w2, 104, 150, 7);
% Function: move_search
% Elapsed time is 0.001351 seconds.
% (109, 155) 
% matched-similarity: 58813
% searched-steps: 2

% [i_end1, j_end1, s_end1, S_NaN1, search_steps1, search_path1] = move_search(f_squr, w2, 104, 150, 8);
% Function: move_search
% Elapsed time is 0.000840 seconds.
% (109, 155) 
% matched-similarity: 58813
% searched-steps: 2

% [i_end1, j_end1, s_end1, S_NaN1, search_steps1, search_path1] = move_search(f_squr, w2, 104, 150, 9);
% Function: move_search
% Elapsed time is 0.001989 seconds.
% (109, 155) 
% matched-similarity: 58813
% searched-steps: 2

% [i_end1, j_end1, s_end1, S_NaN1, search_steps1, search_path1] = move_search(f_squr, w2, 104, 150, 10);
% Function: move_search
% Elapsed time is 0.001285 seconds.
% (109, 155) 
% matched-similarity: 58813
% searched-steps: 2

% [i_end1, j_end1, s_end1, S_NaN1, search_steps1, search_path1] = move_search(f_squr, w2, 104, 150, 11);
% Function: move_search
% Elapsed time is 0.001394 seconds.
% (109, 155) 
% matched-similarity: 58813
% searched-steps: 2

% [i_end1, j_end1, s_end1, S_NaN1, search_steps1, search_path1] = move_search(f_squr, w2, 104, 150, 12);
% Function: move_search
% Elapsed time is 0.001034 seconds.
% (109, 155) 
% matched-similarity: 58813
% searched-steps: 2
