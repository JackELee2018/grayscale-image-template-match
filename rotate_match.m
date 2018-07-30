function [target_num, I_match, J_match, theta_match, s_match, match_info, angle_match_accuracy] = rotate_match(im, P, match_paras, w_cell, w_top_cell, steps, threshold_pecentage1, threshold_pecentage2, search_radius1, search_radius2, polarity, multi_mode)
%ROTATE_MATCH Performs rotateing pyramid template-matching.
%   [I_MATCH, J_MATCH, THETA_MATCH, ANGLE_MATCH_ACCURACY] = ROTATE_MATCH(IM, W,
%   P, MATCH_PARAS, W_CELL, STEPS, VARARGIN) performs pyramid rotating
%   template-matching in image, IM, to find target or targets specified by
%   template, W, and gives out the centroid coordinate(s), I_MATCH and J_MATCH,
%   of the matched-targets and rotating angle of targets, THETA_MATCH.
%   ANGLE_MATCH_ACCURACY is accuracy of rotating angle of targets.
%
%   P is layer number of pyramid.
%
%   MATCH_PARAS is a structrue data for passing matching parameters.
%
%   W_CELL is a cell array that stores shrinking templates that have
%   been rotated different angle for rotating template match.
%
%   STEPS is a cell array that stores template rotating steps in each
%   pyramid layer.
% 
%   MULTI_TARGET is a switch of multi-target mode. If there are more than one
%   targets in the image view, the MULTI_TARGET option should be ON.
%   MULTI-TARTET options can choose from belew items, dedault option is
%   enclosed in brace ({}).
% 
%       {'on'}      Multi-target detection function on.
%
%       'off'       Multi-target detection function off. This function is
%                   canceled.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Preparing for rotate_match %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sigle target detection.
% 
% [w, f_squr, f_gray, match_paras, f] = get_template1();
% 
% steps = get_steps(P);
% 
% [w_top_cell, w_cell, match_info] = get_top_templates(w, P, first_step, match_paras);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Multi-target detection.

% f_m = imread('Pic26_VD078SC_RGB_Bar_Ring_Black02_YUY2_15fps.bmp_MODIFED.jpg');
% f2 = rgb2gray(f_m);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start timing.
tic;

% Get seed targets--TOP_MATCH2
[target_array, s_target, target_angle, target_num, im_cell] = top_match2(im, P, w_cell, w_top_cell, threshold_pecentage1, threshold_pecentage2, match_paras, multi_mode);
% [target_array, s_target, target_angle, target_num, im_cell] = top_match2(im, P, w_cell, w_top_cell, threshold_pecentage1, threshold_pecentage2, match_paras, multi_mode, varargin)

s_better = s_target;
% Refine searched locations and angles--ANGLE_MOVE_SEARCH
end_i = {[ ]};
end_j = {[ ]};
end_s = {[ ]};
end_angle = {[ ]};
angles = {[ ]};
direction = {[ ]};
angle_move_steps = {[ ]};
angles_layer = {[ ]};
S_NaN = {[ ]};
search_steps = {[ ]};
search_path = {[ ]};
match_success = {[ ]};

i_next = {[ ]};
j_next = {[ ]};
angle_next = target_angle;

for i_target = 1:target_num
    i_next{i_target} = target_array{i_target}(1);
    j_next{i_target} = target_array{i_target}(2);
end

for layer = 1:(P+1)
    if steps(layer) ~= 0            
        for i_target = 1:target_num
            if isnan(i_next{i_target}) || isnan(j_next{i_target})
                match_success{layer,i_target} = 0;
                end_i{layer}(i_target) = NaN;
                end_j{layer}(i_target) = NaN;
                end_s{layer}(i_target) = NaN;
                end_angle{layer}(i_target) = NaN;
                angles{i_target} = NaN;
                direction{layer, i_target} = NaN;
                angle_move_steps{layer}(i_target) = NaN;
            else
            [match_success{layer, i_target}, end_i{layer}(i_target), end_j{layer}(i_target), end_s{layer}(i_target), end_angle{layer}(i_target), angles{i_target}, direction{layer, i_target}, angle_move_steps{layer}(i_target), match_info]...
                = angle_move_search(im_cell{P+2-layer}, w_cell{P+2-layer}, angle_next{i_target}, steps(layer), i_next{i_target}, j_next{i_target}, search_radius1, match_paras);
%             target_array{i_target}(1,1), target_array{i_target}(1,2)
%             [end_i, end_j, end_s, end_angle, angles, direction, angle_move_steps, match_info] = angle_move_search(im, w, start_angle, angle_step, start_i, start_j, search_radius, match_paras, varargin);

            angles_mat = cell2mat(angles);
            angles_layer{layer,i_target} = angles_mat;
            end
            
            S_NaN{layer} = NaN;
            search_steps{layer}(i_target) = NaN;
            search_path{layer,i_target} = NaN;
        end
        angle_match_accuracy = steps(layer);

    else
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the last detected rotating angle.
% 
%         for i_target = 1:target_num
%             % Get the final refined angles.
%             i_angle = layer - 1;
%             while i_angle > 0 && isnan(end_angle{i_angle}(i_target))
%                 i_angle = i_angle - 1;
%             end
%             if i_angle == 0
%                 angle_next = target_angle{i_target};
% %                 i_final = target_array{i_target}(1,1);
% %                 j_final = target_array{i_target}(1,2);
%             else
%                 angle_next = end_angle{i_angle}(i_target);
% %                 i_final = end_i{i_angle}(i_target);
% %                 j_final = end_j{i_angle}(i_target);
%             end
%         end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            % Debugger code:
%             if isnan(i_final)
%                 i_final = NaN;
%             end
            % Debugger code:
        
        for i_target = 1:target_num
            if isnan(i_next{i_target}) || isnan(j_next{i_target})
                match_success{layer}(i_target) = 0;
                end_i{layer}(i_target) = NaN;
                end_j{layer}(i_target) = NaN;
                end_s{layer}(i_target) = NaN;
                S_NaN{layer} = NaN;
                search_steps{layer}(i_target) = NaN;
                search_path{layer,i_target} = NaN;
            else                
                w_next = imrotate2(w_cell{P+2-layer}, angle_next{i_target}, match_paras);
    %             imgn = imrotate2(f, angle, match_paras, varargin)
                [match_success{layer}(i_target), end_i{layer}(i_target), end_j{layer}(i_target), end_s{layer}(i_target), S_NaN{layer}, search_steps{layer}(i_target), search_path_array]...
                    = move_search(im_cell{P+2-layer}, w_next, i_next{i_target}, j_next{i_target}, search_radius2);
    %             [match_success, i_end, j_end, s_end, S_NaN, search_steps, search_path] = move_search(im, w, i_start, j_start, search_radius, varargin)

                if match_success{layer}(i_target) == 0
                    warning('ROTATE_MATCH:Move_searchOnly',['Move_search failed when P = ', num2str(layer),' and target_num = ',num2str(i_target)'.'])
                end

                % The path coordination will store in a cell of size mxn (m,n > 1).
                % Reshape array, SEARCH_PATH_ARRAY, to matrix of coordination.
                search_path_matrix = (reshape(search_path_array, 2, length(search_path_array)/2))';
                search_path{layer,i_target} = search_path_matrix;

                match_info.move_search_step = search_steps;
                match_info.move_search_path = search_path;
            end

            end_angle{layer}(i_target) = NaN;
            angles{layer}(1,i_target) = NaN;
            direction{layer}(i_target) = NaN;
            angle_move_steps{layer}(i_target) = NaN;
        end
    end
    
    % Calculate seed-points and angles of next pyramid layer.
    for i_target = 1:target_num
        i_next{i_target} = 2*end_i{layer}(i_target);
        j_next{i_target} = 2*end_j{layer}(i_target);  
        
%     for i_target = 1:target_num
%         if isnan(end_i{layer}(i_target))
%             i_next{i_target} = (2^layer)*target_array{i_target}(1);
%             j_next{i_target} = (2^layer)*target_array{i_target}(2);
%         else
%             i_next{i_target} = 2*end_i{layer}(i_target);
%             j_next{i_target} = 2*end_j{layer}(i_target);
%         end
%         angle_next{i_target} = end_angle{layer}(i_target);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    % Find the last detected rotating angle.

    for i_target = 1:target_num
        % Get the final refined angles.
        i_angle = layer;
        while i_angle > 0 && isnan(end_angle{i_angle}(i_target))
            i_angle = i_angle - 1;
        end
        if i_angle == 0
            angle_next{i_target} = target_angle{i_target};
%                 i_final = target_array{i_target}(1,1);
%                 j_final = target_array{i_target}(1,2);
        else
            angle_next{i_target} = end_angle{i_angle}(i_target);
%                 i_final = end_i{i_angle}(i_target);
%                 j_final = end_j{i_angle}(i_target);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Matching similarity value degenration warning.
    switch polarity
        case 'minimum'
            sign = 1;
        case 'maximim'
            sign = -1;
        otherwise
            error('Unknown type of polarity--rotate_match(109).')
    end
    
    [w_m, w_n] = size(w_cell{layer});
    if layer == P+1
        [w_m0, w_n0] = size(w_cell{layer});
    else
        [w_m0, w_n0] = size(w_cell{layer+1});
    end
    
    for i_target = 1:target_num
        if (end_s{layer}(i_target)/(w_m*w_n)*sign) < (s_better{i_target}/(w_m0*w_n0)*sign)
            warning('ROTATE_MATCH:degenerationWarning', ['Similarity value degenerate in P = ',num2str(layer),' and target = ',num2str(i_target),'.'])
        end
    end
end

% Output.
I_match = end_i;
J_match = end_j;
theta_match = end_angle;
s_match = end_s;

% Results display.
disp('Function: rotate_match running results:');
disp(['Matched targets number: ',num2str(target_num)]);
for layer = 1:(P+1)
    disp(['Image pyramid layer No. ', num2str(layer)]);
    for i_target = 1:target_num
        disp(['    Target No. ', num2str(i_target)]);
        disp(['    Matching coordinate: (', num2str(I_match{layer}(i_target)), ',', num2str(J_match{layer}(i_target)),')']);
        disp(['    Matching angle: ', num2str(theta_match{layer}(i_target))]);
        disp(['    Matching similarity: ', num2str(s_match{layer}(i_target))]);
    end
end
disp(['The final angle match accuracy is: ', num2str(angle_match_accuracy)]);

disp('Function: Rotae_match--total consuming time');
% End timing.
toc;

% Mark targets.
f_marked = im;

for i_target = 1:target_num
    for mark_width = 1:3
        f_marked = mark_target(f_marked, w_cell{1}, i_next{i_target}/2, j_next{i_target}/2, angle_next{i_target}, mark_width, 'red');
%         f_marked = mark_target(f_marked, w_cell{1}, I_match{5}(i_target), J_match{5}(i_target),  theta_match{5}(i_target), mark_width, 'red');
    end
end

figure, imshow(f_marked, [ ]);

for i_target = 1:target_num
    text(j_next{i_target}/2, i_next{i_target}/2,'o','color','r');
    text(j_next{i_target}/2+10, i_next{i_target}/2+10,['Target(',num2str(i_target),'):(',num2str(I_match{5}(i_target)),',',num2str(J_match{5}(i_target)),')'],'color','b');
%     text(J_match{5}(i_target), I_match{5}(i_target),'o','color','r');
%     text(J_match{5}(i_target)+10, I_match{5}(i_target)+10,['Target(',num2str(i_target),'):(',num2str(I_match{5}(i_target)),',',num2str(J_match{5}(i_target)),')'],'color','b');
end


% target_num
% target_array
% target_angle
% s_target
% 
% target_num
% I_match, J_match, theta_match, s_match, match_info, angle_match_accuracy





