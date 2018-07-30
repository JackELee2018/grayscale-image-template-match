function r = GaussianFilter3(im, varargin)
%GAUSSIANFILTER3 Performs the Gaussian Filtering using function, IMFILTER, in MATLAB
%Image Processing Toolbox as a contrastive function to GARSSIANFILTER2
%about time consuming.

%   R = GAUSSIANFILTER3(IM, DELTA, RADIUS, 'replicate') filters
%   the image,IM, with
%   Gaussian Approximation filter whose radius is RADIUS and standard
%   deviation is DELTA. During filtering, the size will be extended by
%   copying intensity values in the boundary in 'replicate' mode and by
%   mirroring the image against the boundary in the 'symmetric' mode.
%
%   R = GAUSSIANFILTER3(IM) filter the image, IM, with default values, DELTA
%   = 1 and RADIUS = 5, using 'replicate' mode.

%   R = GAUSSIANFILTER3(IM, DELTA, RADIUS, 'FIXED', B) pads the image with
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
    
% Let GaussianFilter3 have the same output with GaussianFilter.
    im=double(im);
    
% Pad the input image with specified mode.
switch mode
    case 'replicate'
        r = imfilter(im, GaussianSmooth, 'conv', 'replicate');
    case 'fixed'
        r = imfilter(im, GaussianSmooth, 'conv', b);
    case 'symmetric'
        r = imfilter(im, GaussianSmooth, 'conv', 'symmetric');
    otherwise
        error('Wrong padding option has been inputted.');
end

% End timing.
toc;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the contrast method is not fair. Gaussianfilter3 has packed the same
% function two times, it has more extra time consuming.

% r3 = GaussianFilter2(f);
%%%%%%% Elapsed time is 0.194094 seconds.%%%%%%%

% r4 = GaussianFilter3(f);
%%%%%%% Elapsed time is 0.564885 seconds.%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Without package of extra functions. The IPT function performs more
% excellently. Really nice! There should be imperovements in my programs.

% tic
% w2 = GaussianMatrix(1, 5);
% radius = 5;
% delta = 1;
% im=double(f);
% [m,n]=size(f);
% impad = zeros(m + 2*radius, n + 2*radius);
% r=zeros(m,n);
% impad((radius + 1):(radius + m), (radius + 1):(radius + n)) = im;
% impad((radius + 1):(radius + m), 1:radius) = repmat(im(:, 1), 1, radius);
% impad((radius + 1):(radius + m), (radius + n + 1):(2*radius + n)) = repmat(im(:, n), 1, radius);
% impad(1:radius, :) = repmat(impad(radius + 1, :), radius, 1);
% impad((radius + m + 1):(2*radius + m), :) = repmat(impad((radius + m), :), radius, 1);
% for i = 1:m
%     for j = 1:n
%         r(i, j) = sum(sum(impad(i:(i + 2*radius), j:(j + 2*radius))...
%             .*w2));
%     end
% end
% toc;
%%%%%%% Elapsed time is 0.212610 seconds.%%%%%%%


% tic
% w = GaussianMatrix(1, 5);
% r5 = imfilter(f, w, 'conv', 'replicate');
% toc;
%%%%%%% Elapsed time is 0.065452 seconds.%%%%%%%