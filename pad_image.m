function [im_padded, padding_h, padding_w] = pad_image(im, w, varargin)
%PAD_IMAGE pads image, IM, according to size of template, W. PAD_IMAGE is
%used before carry out template-matching in image, IM. The function
%transforms from function CALCU_SIMILARITY.
%
%   IM_PADDED = PAD_IMAGE(IM, W, VARARGIN). The input image should be gray.
%   Template, W, should be smaller then input image.
%
%   IM_PADDED = PAD_IMAGE(IM, W, PADDING_MATHOD, B).
%   PADDING_MATHOD specifies way of padding the input 
%   image's boundaries. PADDING_METHOD is a string that can have one of the
%   following values. The default value is enclosed in braces ({}). The
%   padding size at each boundary is the half size of the template.
% 
%       {'replicate'}   Pad the boundary by copying pixels at the boundary.
% 
%       'symmetric'     Extend image by symmetrically mirroring boundary image
%                       over the boundary edges.
% 
%       'fixed'         Pad image with a constant, B, that following
%                       PADDING_METHOD. If no B inputed, default value for
%                       B is 0.
% 
%       'im_mean'       Pad image with input image's mean intensity,
%                       im_mean. Its function is the same with using
%                       'fixed' option with B = im_mean.
% 
%       'w_mean'        Pad image with template image's mean intensity.

% Start timing.
tic;
% Verify the right number of inputs.
error(nargchk(2, 4, nargin));

% Store the format of input image for use later.
classin1 = class(im);
classin2 = class(w);

% Transform class of input image to 'doublt' that MATLAB need for
% calculation.
im = double(im);
w = double(w);

% Set default values.
padding_method = 'replicate';
b = 0;

% Predefine for large variables and get input values.
if nargin > 3
    b = varargin{2};
end
if nargin > 2
    padding_method = varargin{1};
end


% Deal with special inputs.
[M, N, imd3] = size(im);
[m, n, wd3] = size(w);
if imd3 ~= 1 || wd3 ~= 1
    error('The program only process gray image for now.');
end
if m < 1 || n < 1
    error('The template rows and columns should not less than 1.');
end
if M <= 2 || N <= 2
    error('The target image should larger than size of 2-by-2.');
end

% Define variables.
NUM = M*N;
num = m*n;

padding_h = floor(m/2);
padding_w = floor(n/2);

im_padded = 255*ones(M + 2*padding_h, N + 2*padding_w);
im_padded(padding_h+1:padding_h+M, padding_w+1:padding_w+N) = im;

im_mean = sum(im(:))/NUM;
w_mean = sum(w(:))/num;

%%%%%%%%%%%%%%%%%%%%%% Padding Image %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pad input image, IM, according to size of template, W.
switch padding_method
    case 'replicate'
        im_padded(1:padding_h, padding_w+1:padding_w+N) = ...
            repmat(im_padded(padding_h+1, padding_w+1:padding_w+N), padding_h, 1);
        im_padded(1:padding_h+M, padding_w+N+1:2*padding_w+N) = ...
            repmat(im_padded(1:padding_h+M, padding_w+N), 1, padding_w);
        im_padded(padding_h+M+1:2*padding_h+M, padding_w+1:2*padding_w+N) = ...
            repmat(im_padded(padding_h+M, padding_w+1:2*padding_w+N), padding_h, 1);
        im_padded(1:2*padding_h+M, 1:padding_w) = ....
            repmat(im_padded(1:2*padding_h+M, padding_w+1), 1, padding_w);
        
    case 'symmetric'
        im_padded((padding_h + 1):(padding_h + M), 1:padding_w) = im(:, (padding_w + 1):-1:2);
        im_padded((padding_h + 1):(padding_h + M), (padding_w + N + 1):(2*padding_w + N)) = im(:, (N - 1):-1:(N - padding_w));
        im_padded(1:padding_h, :) = im_padded((2*padding_h + 1):-1:(padding_h + 2), :);
        im_padded((padding_h + M + 1):(2*padding_h + M), :) = im_padded((padding_h + M - 1):-1:M, :);        
%         error('Sorry! The symmetric padding function is not finished.');
        
    case 'fixed'
        im_padded(1:padding_h, padding_w+1:padding_w+N) = b;
        im_padded(1:padding_h+M, padding_w+N+1:2*padding_w+N) = b;
        im_padded(padding_h+M+1:2*padding_h+M, padding_w+1:2*padding_w+N) = b;
        im_padded(1:2*padding_h+M, 1:padding_w) = b;
        
    case 'im_mean'
        im_padded(1:padding_h, padding_w+1:padding_w+N) = im_mean;
        im_padded(1:padding_h+M, padding_w+N+1:2*padding_w+N) = im_mean;
        im_padded(padding_h+M+1:2*padding_h+M, padding_w+1:2*padding_w+N) = im_mean;
        im_padded(1:2*padding_h+M, 1:padding_w) = im_mean;
        
    case 'w_mean'
        im_padded(1:padding_h, padding_w+1:padding_w+N) = w_mean;
        im_padded(1:padding_h+M, padding_w+N+1:2*padding_w+N) = w_mean;
        im_padded(padding_h+M+1:2*padding_h+M, padding_w+1:2*padding_w+N) = w_mean;
        im_padded(1:2*padding_h+M, 1:padding_w) = w_mean;
        
    otherwise
        error('Unknown padding method has been inputed.(PAD_IMAGE)');
end

disp('Function: pad_image; total elapsed time.');
% End timing.
toc;



        