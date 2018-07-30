function [match_success, end_i, end_j, end_s, end_angle, angles, direction, angle_move_steps, match_info] = angle_move_search(im, w, start_angle, angle_step, start_i, start_j, search_radius, match_paras, varargin)
%ANGLE_MOVE_SEARCH Performs rotating template matching in image, IM, to
%find the best matching angle and point of target specified by template,
%W. Angle search starts from given angle, START_ANGLE. ANGLE_MOVE_SEARCH
%use function MOVE_SEARCH to find best matching point and optimal
%similarity when in a given angle.
% 
%   VARARGIN{1} = POLARITY stands for polarity of similarity metric and can
%   be chosen from below, the default is enclosed in brace ({}).
%       
%       {'minimum'}     The smaller the similarty metric, the better the
%                       matching.
%       'maximum'       The bigger the similarity metric, the better the
%                       matching.
% 
% 
% 
% 
% 
% 
%   Example: [end_i, end_j, end_s, end_angle, angles, direction, angle_move_steps, match_info] = angle_move_search(im, w, start_angle, angle_step, start_i, start_j, search_radius, match_paras, varargin)

% Start timing.
tic;

% Check the right number of inputs.
error(nargchk(8, 9, nargin));

% Set default inputs.
polarity = 'minimum';

if nargin > 8
    polarity = varargin{1};
end

% Judge the angle_search deriction using three test angle.
w_start = imrotate2(w, start_angle, match_paras);
w_forward = imrotate2(w, start_angle + angle_step, match_paras);
w_backward = imrotate2(w, start_angle - angle_step, match_paras);

[match_success1, i1(1), j1(1), s1(1), S1{1}, search_step1(1), search_path1{1}] = move_search(im, w_start, start_i, start_j, search_radius);
[match_success2, i1(2), j1(2), s1(2), S1{2}, search_step1(2), search_path1{2}] = move_search(im, w_forward, start_i, start_j, search_radius);
[match_success3, i1(3), j1(3), s1(3), S1{3}, search_step1(3), search_path1{3}] = move_search(im, w_backward, start_i, start_j, search_radius);

if match_success1 == 0 || match_success2 == 0 || match_success3 == 0
    warning('ANGLE_MOVE_SEARCH:JudgeDeriction',['Move_search failed when angle_move_search when start_angle = ',num2str(start_angle),' angle_step = ',num2str(angle_step),' and start_i = ',num2str(start_i),' and start_j = ',num2str(start_j),'.']);
    match_success = 0;
    end_i = NaN;
    end_j = NaN;
    end_s = NaN;
    end_angle = NaN;
    angles = NaN;
    direction = NaN;
    angle_move_steps = NaN;
    match_info = NaN;
    return
end

match_success = {match_success1};
i = {i1(1)};
j = {j1(1)};
s = {s1(1)};
% S = {S1{1}};
S = S1(1);
% S = S1{1};
search_steps = {search_step1(1)};
% search_path = {search_path1{1}};
search_path = search_path1(1);

angles = {start_angle};

switch polarity
    case 'minimum'
        sign = 1;
    case 'maximum'
        sign = -1;
    otherwise
        error('Unknown polarity option.(angle_move_search)');
end

if s1(1) > s1(2) && s1(1) > s1(3)
    situation1 = 1;
elseif s1(1) < s1(2) && s1(1) < s1(3)
    situation1 = -1;
elseif s1(2) < s1(1) && s1(1) < s1(3)
    situation1 = 2;
elseif s1(2) > s1(1) && s1(1) > s1(3)
    situation1 = -2;
else
    error('Function: ANGLE_MOVE_SEARCH--The angle_step may be too small and invalid.');
end

result = sign * situation1;

switch result
    case 1
        warning('ANGLE_MOVE_SEARCH:AngleDirectionDetect',['Move_search failed when angle_move_search when start_angle = ',num2str(start_angle),' angle_step = ',num2str(angle_step),' and start_i = ',num2str(start_i),' and start_j = ',num2str(start_j),'.']);
        error('Function: ANGLE_MOVE_SEARCH failed.--The angle search failed due to rotating forward and rotating backward are both better than the present angle. Can not find a proper rotating deriction. The angle_move_step may be too small.');
    case -1
        end_i = i1(1);
        end_j = j1(1);
        end_s = s1(1);
        end_angle = start_angle;
        direction = NaN; 
        angle_move_steps = 0;
        
        match_info.move_search_steps = search_step1;
        match_info.move_search_path = search_path1;
        
        % Convert cell to array.
        match_success = cell2mat(match_success);
        angles = cell2mat(angles);
%         end_angle = cell2mat(end_angle);
        
        match_info.S = NaN;
        match_info.move_search_steps = NaN;
        match_info.move_search_path = NaN;
        
        return
        
    case 2
        direction = 'anticlockwise';
        
        match_success{2} = match_success2;
        i{2} = i1(2);
        j{2} = j1(2);
        s{2} = s1(2);
        S{2} = S1{2};
        search_steps{2} = search_step1(2);
        search_path{2} = search_path1{2};
        
        angles{2} = start_angle + angle_step;
        
    case -2
        direction = 'clockwise';
        angle_step = -angle_step;
        
        match_success{2} = match_success3;
        i{2} = i1(3);
        j{2} = j1(3);
        s{2} = s1(3);
        S{2} = S1{3};
        search_steps{2} = search_step1(3);
        search_path{2} = search_path1{3};
                
        angles{2} = start_angle + angle_step;
        
    otherwise
        error('Error case of function ANGLE_MOVE_SEARCH.');
end

% Start angle-move searching.
found = 1;
t = 3;

while(found)
    % Next angle.
    angles{t} = angles{t-1} + angle_step;
    
    % Get new template.
    w_r = imrotate2(w, angles{t}, match_paras);
    [match_success{t}, i{t}, j{t}, s{t}, S{t}, search_steps{t}, search_path{t}] = move_search(im, w_r, start_i, start_j, search_radius);
    
    if match_success{t} == 0
        error(['Move_search failed when (i,j)= ', num2str(start_i),',',num2str(start_j),') and template angle = ',num2str(angles{t})'.'])
    end
    
    if s{t} < s{t-1}
        situation2 = 1;
    elseif s{t} > s{t-1}
        situation2 = -1;
    else
        error('Function: ANGLE_MOVE_SEARCH failed.--The angle search failed due to rotating forward and rotating backward are both better than the present angle. Can not find a proper rotating deriction.');
    end
    
    result2 = sign * situation2;
    
    switch result2
        case 1
            t = t + 1;
        case -1
            found = 0;
            
            % Clear data of next time.
            end_i = i{t-1};
            end_j = j{t-1};
            end_s = s{t-1};
            end_angle = angles{t-1};
            
            angles{t} = [ ];
            
            angle_move_steps = t-1;
            s{t} = [ ];
            match_info.S = S;
            search_steps{t} = NaN;
            match_info.move_search_steps = search_steps;
            search_path{t} = NaN;
            match_info.move_search_path = search_path;
            
        otherwise
            error('Error situation of searching--ANGLE_MOVE_SEARCH.');
    end
end

% Output.
% Convert cell to array.
angles = cell2mat(angles);
% end_angle = cell2mat(end_angle);
match_success = cell2mat(match_success);

disp('Function: angle_move_search--total time consuming:');

% End timing.
toc;

% figure, imshow(S1(1), [ ]);
% size(S1{1})
% size(S{1})
            
    



