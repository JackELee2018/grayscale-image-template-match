function [Rp, R_cell, match_paras] = DownSample4(im, varargin)
%DOWNSAMPLE4 Is a data-structure improved version of DOWNSPALE2, which uses
%the cell array to store the size-changing pyramid output images.
%DOWNSAMPLE4 performs image down-sampling specified times and with
%specified method.
%
%   [GD, R_CELL] = DOWNSAMPLE4(IM, VARARGIN) performs half-size down-sampling of image,
%   IM, by specified time using specified algorithm.
%
%   GD is the specified image that down-sampled specified times.
% 
%   R_CELL is a cell array that stores pyramid images of different size.
%   
%   [GD, R_CELL] = DOWNSAMPLE4(IM, P, 'NEAREST', 'TIME_SAVE') re-do downsampling P times
%   and computes gray values in new image by nearest interpolation using a
%   time-save and space-consuming way, against the space-save and
%   time-consuming way when choose 'SPACE_SAVE' option. The last parameter
%   is only working for the 'BILINEAR' option. P should be positive integer.
%   
%   [GD, R_CELL] = DOWNSAMPLE4(IM, P, 'BILINEAR', 'TIME_SAVE')
%   
%   [GD, R_CELL] = DOWNSAMPLE4(IM, P, 'BICUBIC', 'TIME_SAVE')---this function hasn't
%   be finished.
%   
%   DEFAULT: P = 1, 'BILINEAR', 'TIME_SAVE'.
% 
% 
% Example: [Rp, R_cell, match_paras] = DownSample4(im, varargin)


% Start timing.
tic;

% Verify the correct number of input.
error(nargchk(1, 4, nargin));

% Store the class of input for later use.
classin = class(im);

% Change input calss into 'double', which is hope for MATLAB.
im = double(im);

% Set variables.                              
[m, n, iml] = size(im);  

% Even rows and even column will make bilinear interpolation easier.
% So if numbers of row or column are odd, abandon the last row or column.
if round(m/2)~= m/2
    M = m-1;
else
    M = m;
end
if round(n/2) ~= n/2
    N = n-1;
else
    N = n;
end

% Set default values.
p = 1;
intpol_method = 'bilinear';
timeorspace = 'time_save';

if nargin > 3
    timeorspace = varargin{3};
end
if nargin > 2
    intpol_method = varargin{2}; 
end
if nargin > 1
    p = varargin{1};
    if (round(p) ~= p) || (p <= 0)
        error('Parameter, p, should positive integer.');
    end
    % For 'bilinear' mode, the thresholds are ceil(log2(M)) and
    % ceil(log2(N)); For 'nearest' mode the can be round(log2(M)) and
    % round(log2(N)).
    if (p > round(log2(M))) || (p > round(log2(N)))
        error(['Parameter, p, is too big. P should be at the range of [round(log2(M)), round(log2(N))] = ', '[', num2str(ceil(log2(M))), ', ', num2str(ceil(log2(N))), ']']);
    end
end

   
% Predifine large variables.
R_cell = {im};
% R(M, N) = 125*ones(M, N);     % Pyramid images of the same orginal size without cutoff.
                              % Useful shrinking parts lie on the left-up corner.

                              
% Deal with special inputs.
if iml ~= 1
    error('The program only deals with gray images for now!');
end                              
                              
% Perform the agolithm.
M_new = M;
N_new = N;        
M_old = M;
N_old = N;
match_paras.layer_size_m(1) = M;
match_paras.layer_size_n(1) = N;

switch intpol_method
    case 'nearest'
        for k = 1:p
            M_new = ceil(M_old/2);
            if M_new < 1
                error('The number of new image rows is less than 1! There are too many layers of the pyramid.');
            end
            N_new = ceil(N_old/2); 
            if N_new < 1
                error('The number of new image columns is less than 1! There are too many layers of the pyramid.');
            end              
            R_cell{k+1} = R_cell{k}(1:2:M_old, 1:2:N_old);
            M_old = M_new;
            N_old = N_new;
            match_paras.layer_size_m(k+1) = M_new;
            match_paras.layer_size_n(k+1) = N_new;
        end

    case 'bilinear'      
        switch timeorspace
            case 'time_save'

                    im_cell1 = {zeros(M, N)};
                    im_cell2 = {zeros(M, N)};
                    im_cell3 = {zeros(M, N)};
                    im_cell4 = {zeros(M, N)};
                
                for k = 1:p
                    % Error!!! M and N should be largest even rows and 
                    % column of last new created image.
                    % Processing M and N after images processing.
                   
                    % Transfer the four-neighborhood pixels into four
                    % sub-image.  
                    
                    % Old and new down-sampling image size.
                    M_new = floor(M_old/2);
                    if M_new < 1
                        error('The number of new image rows is less than 1! There are too many layers of the pyramid.');
                    end                    
                    N_new = floor(N_old/2);
                    if N_new < 1
                        error('The number of new image columns is less than 1! There are too many layers of the pyramid.');
                    end                      
                    
                    im_cell1{k+1} = R_cell{k}(1:2:M_old, 1:2:N_old);
                    im_cell2{k+1} = R_cell{k}(2:2:M_old, 1:2:N_old);
                    im_cell3{k+1} = R_cell{k}(1:2:M_old, 2:2:N_old);
                    im_cell4{k+1} = R_cell{k}(2:2:M_old, 2:2:N_old);
                    
                    % Bilinear interpolation.
                    R_cell{k+1} = 0.25*(im_cell1{k+1} + im_cell2{k+1} +...
                        im_cell3{k+1} + im_cell4{k+1});
                    
                    % Shrink M and N.
                    M_old = M_new;
                    N_old = N_new;
                    match_paras.layer_size_m(k+1) = M_new;
                    match_paras.layer_size_n(k+1) = N_new; 
                    
                    % Check odevity of M_old and N_old.
                    if round(M_old/2)~= M_old/2
                        M_old = M_old-1;
                    end
                    if round(N_old/2) ~= N_old/2
                        N_old = N_old-1;
                    end                     
                end
                
            case 'space_save'
                % Predefine first.
                
                for k = 1:p
                    % Check odevity of M_old and N_old. No use for the
                    % first running.
                    if round(M_new/2)~= M_new/2
                        M_new = M_new-1;
                    end
                    if round(N_new/2) ~= N_new/2
                        N_new = N_new-1;
                    end                      
                    % Shrink M and N.
                    M_old = M_new;
                    N_old = N_new;                    
                    for i = 1:2:M_old
                        for j = 1:2:N_old
                            R_cell{k+1}((0.5*(i - 1) + 1), (0.5*(j - 1) + 1)) =...
                                0.25*(R_cell{k}(i, j) + R_cell{k}(i + 1, j) +...
                                R_cell{k}(i, j + 1) + R_cell{k}(i + 1, j +1));                             
                        end
                    end
                    
                    % Shrink M and N.
                    M_new = floor(M_old/2);
                    if M_new < 1
                        error('The number of new image rows is less than 1! There are too many layers of the pyramid.');
                    end                      
                    N_new = floor(N_old/2);
                    if N_new < 1
                        error('The number of new image columns is less than 1! There are too many layers of the pyramid.');
                    end  
                    match_paras.layer_size_m(k+1) = M_new;
                    match_paras.layer_size_n(k+1) = N_new;                    
                end
        end
end

% Convert class and output.
% Rc = R(1:M_new, 1:N_new, p + 1);
Rp = R_cell{p+1};
Rp = im2uint8(mat2gray(Rp));

% figure, imshow(Rp, [ ]);
% for k = 1:p + 1
% figure, imshow(R(:, :, k), [ ]);
% end

disp(['Function: DownSample4; with "', timeorspace, '" mode.']);

% End timing.
toc;

disp(['Size of input image is: ', num2str(size(im))]);
disp(['Size of output image is: ', num2str(size(Rp))]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% subplot %%%%%%%%%%%%%%%%%%%%%%%%%%

% figure;
% left = 0;
% bottom = 0;
% [M, N] = size(R_cell{1});
% % sum_width = 0;
% for plot = 1:5
%     
%     [m, n] = size(R_cell{plot});
%     width = 0.27*(1/N)*n;
% %     sum_width = sum_width + width;
%     height = 1;
% %     height = 0.27*(1/M)*m;
%     subplot(1,5, plot, 'position', [left, bottom, width, height]);
% %     subplot(1,5, plot);
%     left = left + width + 0.05;
%     imshow(R_cell{plot}, [ ]);
%     hold on;
% %     xlim([0, 1024]);
% %     ylim([0, 768]);
% 
% 
% % hold on;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% m = 0.1;
% for plot = 1:5
% subplot(1,5, plot, 'Position', [m,.5, .3, .3]);
% imshow(R_cell{plot}, [ ]);
% % [m, ~] = size(R_cell{plot});
% xlim([0, 1024]);
% ylim([0, 768]);
% m = m + 0.2/plot;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  figure;
% m = 0;
% for plot = 1:5
% subplot(1,5, plot, 'Position', [m,.0, .25, 1]);
% imshow(R_cell{plot}, [ ]);
% % [m, ~] = size(R_cell{plot});
% xlim([0, 1024]);
% ylim([0, 768]);
% m = m + 0.2/plot;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  figure;
% m = 0;
% for plot = 1:5
% subplot(1,5, plot, 'Position', [m,.0, .25, 1]);
% imshow(R_cell{plot}, [ ]);
% % [m, ~] = size(R_cell{plot});
% xlim([0, 1024]);
% ylim([0, 768]);
% m = m + 0.21/plot;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% m = 0;
% for plot = 1:5
% subplot(1,5, plot, 'Position', [m,.0, .25, 1]);
% imshow(R_cell{plot}, [ ]);
% % [m, ~] = size(R_cell{plot});
% xlim([0, 1024]);
% ylim([0, 768]);
% m = m + 0.19/plot;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure;
% left = [0, .25, .4, .45, .48, .5];
% bottom = [0, 0, 0, 0, 0];
% width = [.4, .15, .10, .06, .04];
% height = [.4, .15, .10, .06, .04];
% % height = [1, 1, 1, 1, 1];
% 
% for plot = 1:5
% subplot(1,5, plot, 'Position', [left(plot), bottom(plot), width(plot), height(plot)]);
% imshow(R_cell{plot}, [ ]);
% xlim([0, 1024]);
% ylim([0, 768]);
% end

