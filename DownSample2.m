function [rd, R, match_paras] = DownSample2(im, varargin)
%DOWNSAMPLE2 Is a modified version of DOWNSPALE2 to clear the odd rows or 
%columns processing error. But not finished.
%DOWNSAMPLE2 performs image down-sampling specified times and with
%specified method.
%
%   [GD, R] = DOWNSAMPLE2(IM, VARARGIN) performs half-size down-sampling of image,
%   IM, by specified time using specified algorithm.
%
%   GD is the specified image that down-sampled specified times.
% 
%   R is array of down-sampling images without cutting the size and with the
%   target images lie in the up-left corner. For example, if you want a
%   image of down-sampling P times; R will contain P images and in the last
%   image contain the useful part that you want. The program cut the (P+1)th
%   image to make output image, GD.
% 
%   R is a m-by-n-by-(P+1) matrix of three dimesion. m-by-n is the size of each images. (P+1)
%   represent the original image and p down-sampling images.
%   
%   [GD, R] = DOWNSAMPLE2(IM, P, 'NEAREST', 'TIME_SAVE') re-do downsampling P times
%   and computes gray values in new image by nearest interpolation using a
%   time-save and space-consuming way, against the space-save and
%   time-consuming way when choose 'SPACE_SAVE' option. The last parameter
%   is only working for the 'BILINEAR' option. P should be positive integer.
%   
%   [GD, R] = DOWNSAMPLE2(IM, P, 'BILINEAR', 'TIME_SAVE')
%   
%   [GD, R] = DOWNSAMPLE2(IM, P, 'BICUBIC', 'TIME_SAVE')---this function hasn't
%   be finished.
%   
%   DEFAULT: P = 1, 'BILINEAR', 'TIME_SAVE'.

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
R(:, :, 1) = im;
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
            R(:, :, (k + 1)) = zeros(m, n);
        end
        for k = 1:p
            M_new = ceil(M_old/2);
            if M_new < 1
                error('The number of new image rows is less than 1! There are too many layers of the pyramid.');
            end
            N_new = ceil(N_old/2); 
            if N_new < 1
                error('The number of new image columns is less than 1! There are too many layers of the pyramid.');
            end              
            R(1:M_new, 1:N_new, k + 1) = R(1:2:M_old,...
                1:2:N_old, k);
            M_old = M_new;
            N_old = N_new;
            match_paras.layer_size_m(k+1) = M_new;
            match_paras.layer_size_n(k+1) = N_new;            
        end

    case 'bilinear'      
        switch timeorspace
            case 'time_save'
                for k = 1:p
                    
                    im1(:, :, (k + 1)) = zeros(M, N);
                    im2(:, :, (k + 1)) = zeros(M, N);
                    im3(:, :, (k + 1)) = zeros(M, N);
                    im4(:, :, (k + 1)) = zeros(M, N);
                    R(:, :, (k + 1)) = zeros(m, n);
                end
                
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
                    
                    im1(1:M_new, 1:N_new, (k + 1)) =...
                        R(1:2:M_old, 1:2:N_old, k);
                    im2(1:M_new, 1:N_new, (k + 1)) =...
                        R(2:2:M_old, 1:2:N_old, k);
                    im3(1:M_new, 1:N_new, (k + 1)) =...
                        R(1:2:M_old, 2:2:N_old, k);
                    im4(1:M_new, 1:N_new, (k + 1)) = ...
                        R(2:2:M_old, 2:2:N_old, k);
                    
                    % Bilinear interpolation.
                    R(1:M_new, 1:N_new, (k + 1)) =...
                        0.25*(im1(1:M_new, 1:N_new, (k + 1)) +...
                        im2(1:M_new, 1:N_new, (k + 1)) +...
                        im3(1:M_new, 1:N_new, (k + 1)) +...
                        im4(1:M_new, 1:N_new, (k + 1)));
                    
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
                    R(:, :, (k + 1)) = zeros(m, n);
                end
                
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
                            R((0.5*(i - 1) + 1), (0.5*(j - 1) + 1), k + 1) =...
                                0.25*(R(i, j, k) + R(i + 1, j, k) +...
                                R(i, j + 1, k) + R(i + 1, j +1, k));                             
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
Rc = R(1:M_new, 1:N_new, p + 1);
rd = im2uint8(mat2gray(Rc));
figure, imshow(rd, [ ]);
% for k = 1:p + 1
% figure, imshow(R(:, :, k), [ ]);
% end

disp(['Function: DownSample4; with "', timeorspace, '" mode.']);

% End timing.
toc;

disp(['Size of input image is: ', num2str(size(im))]);
disp(['Size of output image is: ', num2str(size(rd))]);
