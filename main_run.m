function bSuccess = main_run()
%MAIN_RUN Runs examples of image matching
%
%   BSUCCESS = MAIN_RUN() runs two example matching process.
%   First one is single target matching, second one is multi-targets matching.
%   Reture 1 when function runs successly, otherwise 0.
%

bSuccess = 0;

% FIRST EXAMPLE:
% Get template image
[w, f_squr, f_gray, match_paras, f] = get_template1();

% Set number of layers of image pyramid
match_paras.P = 4;

% Get angle step of rotating templetes
steps = get_steps(4);

% Generate templates using in the top of image pyramid
[w_top_cell, w_cell, match_info] = get_top_templates(w, 4, 45, match_paras);

match_paras.first_step = 45;

% Single target matching
[target_num, I_match, J_match, theta_match, s_match, ...
    match_info, angle_match_accuracy] = rotate_match(f_squr, 4, ...
    match_paras, w_cell, w_top_cell, steps, 0.8, 0.8, 3, 3, 'minimum', 'on')
    
    
% SECOND EXAMPLE-MULTI
f_m = imread('Pic26_VD078SC_RGB_Bar_Ring_Black02_YUY2_15fps.bmp_MODIFED.jpg');
f2 = rgb2gray(f_m);

[target_num, I_match, J_match, theta_match, s_match, match_info, ...
angle_match_accuracy] = rotate_match(f2, 4, match_paras, w_cell, w_top_cell, steps, ...
0.4, 0.4, 6, 6, 'minimum', 'on');