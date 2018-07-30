function f_marked = mark_target(im, w, i_target, j_target, theta_target, frame_offset, varargin)
%MARK_TARGET regions a target in image, IM, using a frame of template size
%that given by template, W, which should be squared.
% 
%	FRAME_OFFSET is offset between marking frame and the boundary of image,
%	IM. Positive frame_offset will get marking frame that bigger than given
%	template.
%   
%	IM can be RGB and gray image.


% Check right number of inputs.
error(nargchk(6, 7, nargin));

% Get inputs and set default values.
color = 'blue';

if nargin > 6
    color = varargin{1};
end

switch color
    case 'blue'
        color_value = [0, 0, 255];
    case 'red'
        color_value = [255, 0, 0];
    case 'green'
        color_value = [0, 255, 0];
    otherwise
        error('Unknown color option has been inputted--Function: mark_target')
end

[m_im, n_im, ~] = size(im);
[m_w, n_w] = size(w);

i_centre = ceil(m_w/2);
j_centre = ceil(n_w/2);

% Get frame's angular points.
top_left_i = i_target - i_centre - frame_offset;
top_left_j = j_target - j_centre - frame_offset;

bottom_left_i = i_target + (m_w - i_centre + 1) + frame_offset;
bottom_left_j = j_target - j_centre - frame_offset;

bottom_right_i = i_target + (m_w - i_centre + 1) + frame_offset;
bottom_right_j = j_target + (n_w - j_centre + 1) + frame_offset;

top_right_i = i_target - i_centre - frame_offset;
top_right_j = j_target + (n_w - j_centre + 1) + frame_offset;

% Get squared template frame's coordinates.
frame_coordinates = zeros(2*(2*(frame_offset+1)+m_w) + 2*(2*(frame_offset+1)+n_w) - 4, 2);
frame_coordinates_rotated = zeros(2*(2*(frame_offset+1)+m_w) + 2*(2*(frame_offset+1)+n_w) - 4, 2);

frame_coordinates(1:(2*(frame_offset+1)+m_w), 1) = (top_left_i : bottom_left_i)';  % Octave
frame_coordinates(1:(2*(frame_offset+1)+m_w), 2) = j_target - j_centre - frame_offset;

frame_coordinates((2*(frame_offset+1)+m_w)+1:(2*(frame_offset+1)+m_w)+(2*(frame_offset+1)+n_w)-1, 1) =  i_target + (m_w - i_centre + 1) + frame_offset;
frame_coordinates((2*(frame_offset+1)+m_w)+1:(2*(frame_offset+1)+m_w)+(2*(frame_offset+1)+n_w)-1, 2) = bottom_left_j + 1 : bottom_right_j;

frame_coordinates((2*(frame_offset+1)+m_w)+(2*(frame_offset+1)+n_w):2*(2*(frame_offset+1)+m_w)+(2*(frame_offset+1)+n_w)-2, 1) = bottom_right_i - 1 : -1 : top_right_i;
frame_coordinates((2*(frame_offset+1)+m_w)+(2*(frame_offset+1)+n_w):2*(2*(frame_offset+1)+m_w)+(2*(frame_offset+1)+n_w)-2, 2) = bottom_right_j;

frame_coordinates(2*(2*(frame_offset+1)+m_w)+(2*(frame_offset+1)+n_w)-1:2*(2*(frame_offset+1)+m_w)+2*(2*(frame_offset+1)+n_w)-4, 1) = top_right_i;
frame_coordinates(2*(2*(frame_offset+1)+m_w)+(2*(frame_offset+1)+n_w)-1:2*(2*(frame_offset+1)+m_w)+2*(2*(frame_offset+1)+n_w)-4, 2) = top_right_j - 1 : -1 : top_left_j + 1;

%%%%%%%%%%%%%%%%%%%%%Debugging code%%%%%%%%%%%%%%%%%%%%
% % Display template frame without rotating.
% im_test = im;
% im_test2 = ones(size(im));
% 
% for i = 1:(2*(frame_offset+1)+m_w) + (2*(frame_offset+1)+n_w) - 4
%     if frame_coordinates(i, 1) <= m_im && frame_coordinates(i, 1) > 0 && frame_coordinates(i, 2) <= n_im && frame_coordinates(i, 2) > 0
%         im_test(frame_coordinates(i, 1), frame_coordinates(i, 2)) = 0;
%         im_test2(frame_coordinates(i, 1), frame_coordinates(i, 2)) = 0;
%     end
% end
% 
% figure, imshow(im_test, [ ]);
% figure, imshow(im_test2, [ ]);
%%%%%%%%%%%%%%%%%%%%%Debugging code%%%%%%%%%%%%%%%%%%%%

% Transform rectangular coordinates to polar coordinates.
polar_theta = zeros(1, length(frame_coordinates));
polar_radius = zeros(1, length(frame_coordinates));

for i = 1:2*(2*(frame_offset+1)+m_w) + 2*(2*(frame_offset+1)+n_w) - 4
    x = frame_coordinates(i, 1);
    y = frame_coordinates(i, 2);
    
    % Tangent function is a periodic function with period of PI. When theta
    % is greater than PI radian, add a PI to the calculated theta.
    if x > i_target
        % Pay attation!!! Denominator can't be zere. Add a EPS to every
        % denominator.
        if y > j_target
            
            polar_theta(i) = atan((x - i_target)/(y - j_target + eps));
        else
            
            polar_theta(i) = pi - atan((x - i_target)/(j_target - y + eps));
        end
    else
        if y > j_target
            polar_theta(i) =2*pi - ((atan((x - i_target)/(j_target - y + eps)) + pi) - pi);
%             polar_theta(i) =pi + 2*pi - atan((x - i_target)/(j_target - y + eps)) + pi;
%             polar_theta(i) =4*pi - atan((x - i_target)/(j_target - y + eps));
        else
            
            polar_theta(i) = atan((x - i_target)/(y - j_target + eps)) + pi;
        end
    end
    polar_radius(i) = sqrt((x - i_target)^2 + (y - j_target)^2);
    
end

% Rotate the target frame in polar coordinate.

% Transform degree measure to radian measure.
target_radian = (pi/180)*(-theta_target);

polar_theta_rotated = polar_theta + target_radian;

% Transform rotated polar coordinates back to rectangular coordinates.
for i = 1:(2*(2*(frame_offset+1)+m_w) + 2*(2*(frame_offset+1)+n_w) - 4)
    
    x = ceil(i_target + sin(polar_theta_rotated(i))*polar_radius(i));
    y = ceil(j_target + cos(polar_theta_rotated(i))*polar_radius(i));

    frame_coordinates_rotated(i, 1) = x;
    frame_coordinates_rotated(i, 2) = y;
    
end

% Draw the made frame in the image, IM.

% Make sure the image, IM, is GRB image.
[~, ~, rgb] = size(im);
if rgb == 1
    f_marked(:, :, 1) = uint8(im);
    f_marked(:, :, 2) = uint8(im);
    f_marked(:, :, 3) = uint8(im);
elseif rgb == 3
    f_marked = im;
elseif rgb ~= 3 && rgb ~= 1
    error('Input image should be gray or RGB image--Function:mark_target.')
end

for i = 1:2*(2*(frame_offset+1)+m_w) + 2*(2*(frame_offset+1)+n_w) - 4
    if frame_coordinates_rotated(i, 1) <= m_im && frame_coordinates_rotated(i, 1) > 0 && frame_coordinates_rotated(i, 2) <= n_im && frame_coordinates_rotated(i, 2) > 0
        f_marked(frame_coordinates_rotated(i, 1), frame_coordinates_rotated(i, 2),:) = color_value;
    end
end

% Output.
% figure, imshow(f_marked, [ ]);
        




