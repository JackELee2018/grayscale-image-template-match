function r = GaussianFilter2(im, varargin)
%GAUSSIANFILTER Performs Gaussian Approximation filtering for Image Pryamid
%algorithm using a square filter that create by function GaussianMatrix.
%The image should be gray. It does the same job with function
%GaussianFilter,but in a time-saving way. The algorithm here will extend
%and pad the original image before carry out filtering.
%This code is far more efficient than the function GAUSSIANFILTER, a
%Internet code, but less efficient than the IPT function IMFILTER which use
%a embedded C code to execute the filtering.
%GAUSSIANFILTER3 is countpart function about time-consuming, which use the 
%former IPT function, IMFILTER.
%
%   R = GAUSSIANFILTER2(IM, DELTA, RADIUS, 'replicate') filters
%   the image,IM, with
%   Gaussian Approximation filter whose radius is RADIUS and standard
%   deviation is DELTA. During filtering, the image size will be extended by
%   copying intensity values in the boundary when choose 'replicate' mode. 
%   And when choose the 'symmetric' mode, it mirror the boundary intensity
%   to extend the image.
%
%   R = GAUSSIANFILTER2(IM) filter the image, IM, with default values, DELTA
%   = 1 and RADIUS = 5, using 'replicate' mode.

%   R = GAUSSIANFILTER2(IM, DELTA, RADIUS, 'FIXED', B) pads the image with
%   fixed value, B. The default value for B is 0.

% Start timing.
tic

% Check the number of inputs.
error(nargchk(1, 5, nargin));

% Set default values.
b = 0;
mode = 'replicate';
radius = 5;
delta = 1;

if nargin > 4
    b = varargin{4};
end
if nargin > 3
    mode = varargin{3};
end
if nargin > 2
    radius = varargin{2};
end
if nargin > 1
    delta = varargin{1};
end

% Get the Gaussian Approximation Matrix.
    GaussianSmooth=GaussianMatrix(delta,radius);
    
% For making output and input have the same data format.
    classin = class(im);
    
    im=double(im);
    [m,n]=size(im);
% IMPAD stores the padded image.
    impad = zeros(m + 2*radius, n + 2*radius);
    r=zeros(m,n);
    
% Pad the input image with specified mode.
switch mode
    case 'replicate'
% Create the padded image, IMAPD, with specified mode.

% Error assignment: 'impad((radius + 1):(radius + m), 1:radius) = im(:, 1)';
% left side of equal sign is matrix of size (m-by-(radius + m)); then right
% side of equal sign is a array of size (m-by-1). Use the function REPMAT.
        impad((radius + 1):(radius + m), (radius + 1):(radius + n)) = im;
        impad((radius + 1):(radius + m), 1:radius) = repmat(im(:, 1), 1, radius);
        impad((radius + 1):(radius + m), (radius + n + 1):(2*radius + n)) = repmat(im(:, n), 1, radius);
        impad(1:radius, :) = repmat(impad(radius + 1, :), radius, 1);
        impad((radius + m + 1):(2*radius + m), :) = repmat(impad((radius + m), :), radius, 1);

    case 'symmetric'
        impad((radius + 1):(radius + m), (radius + 1):(radius + n)) = im;
        impad((radius + 1):(radius + m), 1:radius) = im(:, (radius + 1):-1:2);
        impad((radius + 1):(radius + m), (radius + n + 1):(2*radius + n)) = im(:, (n - 1):-1:(n - radius));
        impad(1:radius, :) = impad((2*radius + 1):-1:(radius + 2), :);
        impad((radius + m + 1):(2*radius + m), :) = impad((radius + m - 1):-1:m, :);
    case 'fixed'
        impad((radius + 1):(radius + m), (radius + 1):(radius + n)) = im;
        impad((radius + 1):(radius + m), 1:radius) = b;
        impad((radius + 1):(radius + m), (radius + n + 1):(2*radius + n)) = b;
        impad(1:radius, :) = b;
        impad((radius + m + 1):(2*radius + m), :) = b;
    otherwise
        error('Wrong padding mode has been inputted.');
end
        
% perform the filtering with specified mode.
for i = 1:m
    for j = 1:n
        r(i, j) = sum(sum(impad(i:(i + 2*radius), j:(j + 2*radius))...
            .*GaussianSmooth));
    end
end

% Convert R to the class of the input image.
%r = changeclass(classin, r);

% End timing.
toc;

    
    
    