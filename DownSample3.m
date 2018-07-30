function rd = DownSample4(im, varargin)
%DOWNSAMPLE2 Performs pixels down-sampling with specified method.
%
%   GD = DOWNSAMPLE2(IM, VARARGIN) performs half-size down-sampling of image,
%   IM, using specified algorithm.
%   
%   GD = DOWNSAMPLE2(IM, 'NEAREST', 'TIME_SAVE') downsamples the image
%   and computes gray values in new image by nearest interpolation using a
%   time-save and space-consuming way, against the space-save and
%   time-consuming way when choose 'SPACE_SAVE' option. The last parameter
%   work only when 'NEAREST' option.
%   
%   Pay attation! Size of GD may be different when use 'NEAREST' option and
%   'BILINEAR' option. Size of GD in 'BILINEAR' mode will be one row or one
%   column less than that in 'NEAREST'.
% 
%   GD = DOWNSAMPLE2(IM, 'BILINEAR', 'TIME_SAVE')
%   
%   GD = DOWNSAMPLE2(IM, 'BICUBIC', 'TIME_SAVE')---this function hasn't
%   be finished.
%   
%   DEFAULT: 'BILINEAR', 'TIME_SAVE'.

% Start timing.
tic;

% Verify the correct number of input.
error(nargchk(1, 3, nargin));

% Store the class of input for later use.
classin = class(im);

% Change input calss into 'double', which is hope for MATLAB.
im = double(im);

% Set default values.
intpol_method = 'bilinear';
timeorspace = 'time_save';

if nargin > 2
    timeorspace = varargin{2};
end
if nargin > 1
    intpol_method = varargin{1}; 
end

% Predifine large variables.

% Define variables.
[m, n, iml] = size(im);

% Deal with special inputs.
if iml ~= 1
    error('The program only deals with gray images for now!');
end

% Even rows and even column will make bilinear interpolation easier.
% So if numbers of row or column are odd, abandon the last row or column.
if round(m/2)~= m/2
    m_even = m-1;
else
    m_even = m;
end
if round(n/2) ~= n/2
    n_even = n-1;
else
    n_even = n;
end


% Perform the agolithm.algorithm.
switch intpol_method
    case 'nearest'
%         tic;
        rd = im(1:2:end, 1:2:end);
%         toc;
        
    case 'bilinear'
        switch timeorspace
            case 'time_save'
%                 tic;
                % Transfer the four-neighborhood pixels into four
                % sub-image.
                % Error!!! When m or n are old numbers, im1, im2, im3 and
                % im4 will has different size. The following calculation
                % fail! So change m and n to even numbers. But, this will
                % make size of output images of 'bilinear' option one row
                % or one column less than that of 'nearest' option.
%                 im1 = im(1:2:end, 1:2:end);
%                 im2 = im(2:2:end, 1:2:end);
%                 im3 = im(1:2:end, 2:2:end);
%                 im4 = im(2:2:end, 2:2:end);
                
                im1 = im(1:2:m_even, 1:2:n_even);
                im2 = im(2:2:m_even, 1:2:n_even);
                im3 = im(1:2:m_even, 2:2:n_even);
                im4 = im(2:2:m_even, 2:2:n_even);
% Debugging code.                
%                 whos im1 im2 im3 im4
%                 minmax2(im1);
%                 minmax2(im2);
%                 minmax2(im3);
%                 minmax2(im4);
%                 figure, imshow(im1, [0, 255]);
%                 figure, imshow(im2, [0, 255]);
%                 figure, imshow(im3, [0, 255]);
%                 figure, imshow(im4, [0, 255]);

                % Bilinear interpolation.
                rd = 0.25*(im1 + im2 + im3 + im4);
%                 toc;

            case 'space_save'
%                 tic;
                
                for i = 1:2:m_even
                    for j = 1:2:n_even
                        rd((0.5*(i - 1) + 1), (0.5*(j - 1) + 1)) = 0.25*(im(i, j) + ...
                            im(i + 1, j) + im(i, j + 1) + im(i + 1, j +1));
                    end
                end
%                 toc;
                
            otherwise
                error('Wrong Time_Spare mode has been inputed.');
        end
        
    case 'bicubic'
        error('Sorry, the BICUBIC interpolation function has not been done.')
    otherwise
        error('Wrong interpolation mode has been inputed.');
end

% Convert class and output.
rd = im2uint8(mat2gray(rd));
figure, imshow(rd, [ ]);

disp(['Function: DownSample3; with "', timeorspace, '" mode.']);
% End timing.
toc;

disp(['Size of input image is: ', num2str(size(im))]);
disp(['Size of output image is: ', num2str(size(rd))]);

%EXAMPLE%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r1 = DownSample3(f);
% Elapsed time is 0.009187 seconds.
% Elapsed time is 0.130100 seconds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r1 = DownSample3(f, 'bilinear', 'space_save');
% Elapsed time is 0.192923 seconds.
% Elapsed time is 0.303288 seconds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r1 = DownSample3(f, 'nearest');
% Elapsed time is 0.001404 seconds.
% Elapsed time is 0.158995 seconds.