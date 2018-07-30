function [w_top_cell, w_cell, match_info] = get_top_templates(w, P, first_step, match_paras)
%GET_TOP_TEMPLATES Creates rotating templates with equal rotating angle
%difference, which used in top layer pyramid-rotating-matching processing.
%GET_TOP_CELL usually runs before function, ROTATE_MATCH to prepare datum
%for ROTATE_MATCH, which will save time of matching processing.
% 
%   FIRST_STEP is rotating angle step of templates in W_TOP_CELL.
%   FIRST_STEP shouldn't be too small to increase consuming time, neither
%   too big to miss right targets. Refer to "Seventh Experiment_rotation
%   sensitivity(2)". So far, for right matching, FIRST_STEP should no
%   smaller than 20 degree. And FIRST_STEP should make a circle divided
%   symmetrically to avoid matching 180 degree rotated targets.
%
%   W is original template.
%
%   P is layer number of image pyramid.
%
% 
% Example: [w_top_cell, w_cell, match_paras] = get_top_templates(w, P, first_step, match_paras)

% Clear data may in the variable, W_TOP_CELL.
% w_top_cell = [ ];
% These step should be outside the function and before assigning variable,
% W_TOP_CELL.


% Create template pyramid.
[~, w_cell] = DownSample4(w, P);
template_counter = 1;

% Store origial template and create rotate templates.
w_top_cell = {w_cell{P+1}};
for angle = first_step:first_step:360-first_step
    template_counter = template_counter + 1;
    w_top_cell{template_counter}=imrotate2(w_cell{P+1}, angle, match_paras);
end

% Output.
% close all;
subplot_row = ceil(template_counter/4);
figure;
for plot = 1:template_counter
    subplot(subplot_row, 4, plot); 
%     axis([0 10 0 20]);
    imshow(w_top_cell{plot}, [ ]);
    axis([0 30 0 30]);
end
   
% Transfer paramenters.
match_info.P = P;
match_info.first_step = first_step;


