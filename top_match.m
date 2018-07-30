function [target_array, s_target, target_angle, target_num, im_cell] = top_match(im, P, w_cell, w_top_cell, threshold_pecentage, match_paras, multi_mode, varargin)
%TOP_MATCH performs template matching in the top layer of the image pyramid
%to search all targets in image, IM, and get the rough locations and
%directions of targets.
% 
% 
% 
%   THRESHOLD_PECENTAGE is pecentage of a similarty threshold, which used
%   to judge whether a matched-similarity peak can be a target.
%   THRESHOLD_PECENTAGE should be a decimal and in the form of taking the
%   similarity metric polarity as 'maximum'.
%   
%   MULTI_MODE is a switch of multi-targets detection function. Its value
%   can be one of following items, the default one is enclosed in brace
%   ({}).
%       {'on'}      Multi-targets detection function on.
% 
%       'off'       multi-targets detection function off.
% 
% Example: [target_array, target_angle, target_num, im_cell] = top_match(im, P, w_top_cell, threshold_pecentage, match_paras, multi_mode, varargin) 

% Set default values.
% multi_mode = 'on';
if threshold_pecentage < 0 || threshold_pecentage > 1
    error('Threshold_pecentage should be in range of [0, 1].')
end

% Start timing.
tic;

% Get pyramid of image, IM.
[~, im_cell] = DownSample4(im, P);

% Perform top-matching.
S_top = {[ ]};
I_top = NaN(1, (ceil(360/match_paras.first_step)));
J_top = NaN(1, (ceil(360/match_paras.first_step)));
peak_top = NaN(1, (ceil(360/match_paras.first_step)));

for i_top = 1:length(w_top_cell)
    [S_top{i_top}, I_top(i_top), J_top(i_top), peak_top(i_top), polarity] = calcu_similarity(im_cell{P+1}, w_top_cell{i_top});
end

% Get polarity of similarity metric which CALCU_SIMILARITY assigns to MATCH_PARAS.
% polarity = match_paras.polarity;
% MATCH_PARAS can't be output variable.

switch match_paras.P
    case 1
        threshold_offset = 15000;
    case 2
        threshold_offset = 10000;
    case 3
        threshold_offset = 5000;
    case 4
        threshold_offset = 1000;
    case 5
        threshold_offset = 500;
    case 6
        threshold_offset = 100;
    otherwise
        error('Sorry, other offsets of higher pyramid laryer is not ready--calcu_similarity2.')
end

switch polarity
    case 'minimum'
        s_peak = min(peak_top);
        sign = 1;
        s_threshold = (threshold_offset * (1 - threshold_pecentage)) + s_peak;
    case 'maximum'
        s_peak = max(peak_top);
        sign = -1;
        s_threshold = s_peak - (threshold_offset * (1 - threshold_pecentage));
    otherwise
        error('Unknown polarity of similarity.--top_match');
end

peak_index = find(peak_top == s_peak);

switch multi_mode
    case 'off'
        target_num = 1;
        target_array = ([I_top(peak_index), J_top(peak_index)]);
        s_target = s_peak;
        target_angle = (peak_index - 1) * match_paras.first_step;
        
    case 'on'
        target_num = 1;
        target_array = {[I_top(peak_index), J_top(peak_index)]};
        s_target = {s_peak};
        target_angle = {(peak_index - 1) * match_paras.first_step};
        
        for i_top = 1:length(w_top_cell)
            if (peak_top(i_top)*sign) < (s_threshold*sign)
                target_flag = 1;
                for i_target = 1:target_num
                    % If the distance between new matching point and
                    % existed matching target is larger than width of template
                    % witheout rotating, we think the new matching point is a
                    % new target.
                    if sqrt((target_array{i_target}(1,1)-I_top(i_top))^2 + (target_array{i_target}(1,2)-J_top(i_top))^2) < min(size(w_cell{P+1}))
                        % If distance between new matching point and existed
                        % matching target smaller than width of template,
                        % the two matching points  belongs to the same
                        % target. That comparing their matching similarity,
                        % choose the better.
                        target_flag = 0;
                        % Comparing two matching similarity. Choose better.
                        if (peak_top(i_top)*sign) < (s_target{i_target}*sign)
                            target_array = {[I_top(i_top), J_top(i_top)]};
                            s_target = peak_top(i_top);
                            target_angle = {(i_top - 1) * match_paras.first_step};
                        end
                        % Within range of width of template, there can be 
                        % only one target whose distance between new
                        % matching point is smaller than width of
                        % template. End the loop.
                        break
                    end
                end
                % If the new matching point belongs to no exist points, so
                % add it into the target array.
                if target_flag == 1
                    target_num = target_num + 1;
                    target_array{target_num} = [I_top(i_top), J_top(i_top)];
                    s_target{target_num} = {peak_top(i_top)};
                    target_angle{target_num} = (i_top - 1) * match_paras.first_step;
                end
            end
        end
    otherwise
        error('Unknown option for MULTI_MODE--top_match.')
end

% Output.
disp('Function: top_match--total time consuming:');

% End timing.
toc;

for i_top = 1:length(w_top_cell)
    disp(['Matching point(',num2str(i_top),'): (',num2str(I_top(i_top)),',',num2str(J_top(i_top)),') ',num2str(peak_top(i_top))]);
end



