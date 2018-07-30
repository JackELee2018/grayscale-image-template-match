function [template, f_squr, f_gray, match_paras, f] = get_template1()
%GET_TEMPLATE1 Gets a square chip template from a specified image
%'E:\MatlabPrograms\OriginalChipsImages\Test02_VD087SC_RGB_Monochrome\Pic26_VD078SC_RGB_Bar_Ring_Black02_YUY2_15fps.bmp'
%I have backuped the image as 'motherimage01.bmp' in the same file directory of
%the function code file.
%
%   [TEMPLATE, F_SQUR, F_GRAY, MATCH_PARAS, F] = GET_TEMPLATE1() has no inputs.
%   Example:
%   [w, f_squr, f_gray, match_paras, f] = get_template1();

% tic;

f = imread('Pic26_VD078SC_RGB_Bar_Ring_Black02_YUY2_15fps.bmp');
% figure, imshow(f, [ ]);
f_gray = rgb2gray(f);
figure, imshow(f_gray, [ ]);

% Get background mean intensity.
sub_background = f_gray(21:722, 25:300);
% figure, imshow(sub_background, [ ]);
sub_background_size = size(sub_background);
match_paras.background_mean = sum(sub_background(:))/...
    (sub_background_size(1, 1)*sub_background_size(1, 2));

% match_paras.background_mean =120.2526.

% Rotate image square. Get the rotating angle by picking two chips' points.
theta = atand(21/466);
% >>theta = atand(21/466)
% theta =
%     2.5803
% format long
% theta = atand(21/466)
% theta =
%    2.580252934859366

f_squr = imrotate2(f_gray, theta, match_paras);
disp(['The rotated angle is ', num2str(theta), ' degree.']);
disp(['Size of original image is: ', num2str(size(f))]);
disp(['Size of squared image is: ', num2str(size(f_squr))]);


% Get chip range in image by hand.
template = f_squr(314:563, 365:875);
figure, imshow(template, [ ]);
disp(['Size of template is: ', num2str(size(template))]);

% Output.


% toc;
% disp('Rotate image take main time consuming here.');
