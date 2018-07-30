function r = GaussianFilter(im, varargin)
%GAUSSIANFILTER Performs Gaussian Approximation filtering for
%Image Pyramid algorithm using a square filter that
%creates by function GaussianMatrix.The image can be RGB or gray.
% The code is a modified version of a internet copy:
%http://wenku.baidu.com/link?url=-VI--ZL45KzOwtantrIFoQS3Moi1VNxM57X
%qc6tv9aIUAJ3H1axMw8EXiqET0wwpRIqlGDyy4aFfENBXS7ncXziPpdkHiB1Xo9S-dnF8MXa
%And its efficiency is very low. A running will take more 20 seconds.
%
%   R = GAUSSIANFILTER(IM, DELTA, RADIUS, 'replicate') filters
%   the image,IM, with
%   Gaussian Approximation filter whose radius is RADIUS and standard
%   deviation is DELTA. During filtering, the size will be extended by
%   copying intensity values in the boundary in 'replicate' mode and by
%   mirroring the image against the boundary in the 'symmetric' mode.
%
%   R = GAUSSIANFILTER(IM) filter the image, IM, with default values, DELTA
%   = 1 and RADIUS = 5, using 'replicate' mode.

% Start timing.
tic

% Check the number of inputs.
error(nargchk(1, 4, nargin));

% Set default values.
mode = 'replicate';
delta = 1;
radius = 5;

if nargin < 3
    mode = varargin{3};
end
if nargin < 2
    radius = varargin{2};
end
if nargin < 1
    delta = varargin{1};
end

% Get the Gaussian Approximation Matrix.
    GuassionSmooth=GaussianMatrix(delta,radius);
    
    im=double(im);
    [m,n,z]=size(im);
    r=zeros(m,n,z);
% perform the filtering with specified mode.
switch mode
    case 'replicate'
        for i=1:m
            for j=1:n
                for k=-radius:radius
                    for l=-radius:radius
    % error line:         x1=i+k+1;		change to:
                        x1=i+k;			% Error; Boundary should be ¡®x1 = i + k¡¯.
                                        % --byJackLee, 2016.1.12
                        if x1<1         % Index of Image in MATLAB begins from 1.
    % error line:           x1=x1-2*k+1;		change to:
                            x1 = 1;		% Pixels outside boundary should be 
                                        % ¡®x1 = 1¡¯ corresponding to 
                                        % option ¡®replicate¡¯and¡®x1 = -x1 +2¡¯
                                        % corresponding to ¡®symmetric¡¯. 
                                        % -byJackLee, 2016.1.12
                        end
                        if x1>m				% The last index of Image is m.
    % error line:             x1=x1-2*k-1;		change to:
                            x1= m;          % ¡®Replicate¡¯ should be ¡®x1 = m¡¯;
                                            % ¡®Symmetric¡¯ should be ¡®x1 = x1 - 2¡¯.
                                            % -byJackLee, 2016.1.12
                        end
    % error line:             x2=j+l+1;		change to:
                            x2 = j + l;		% The same as above.
                        if x2<1
    % error line              x2=x2-2*l+1;		change to:
                            x2 = 1;		% The same as above.
                        end
                        if x2>n
    % error line:             x2=x2-2*l-1;		change to:
                            x2 = n;		% The same as about.
                        end
                        r(i,j,:)=r(i,j,:)+...
                            im(x1,x2,:)*GuassionSmooth(k+radius+1,l+radius+1);
                    end
                end
            end
        end
        
    case 'symmetric'
        for i=1:m
            for j=1:n
                for k=-radius:radius
                    for l=-radius:radius
    % error line:         x1=i+k+1;		change to:
                        x1=i+k;			% Error; Boundary should be ¡®x1 = i + k¡¯.
                                        % --byJackLee, 2016.1.12
                        if x1<1         % Index of Image in MATLAB begins from 1.
    % error line:           x1=x1-2*k+1;		change to:
                            x1 = -x1 + 2;		% Pixels outside boundary should be 
                                                % ¡®x1 = 1¡¯ corresponding to 
                                                % option ¡®replicate¡¯and¡®x1 = -x1 +2¡¯
                                                % corresponding to ¡®symmetric¡¯. 
                                                % -byJackLee, 2016.1.12
                        end
                        if x1>m				% The last index of Image is m.
    % error line:             x1=x1-2*k-1;		change to:
                            x1= 2*m - x1;        % ¡®Replicate¡¯ should be ¡®x1 = m¡¯;
                                                 % ¡®Symmetric¡¯ should be ¡®x1 = 2*m - x1¡¯.
                                                 % -byJackLee, 2016.1.12
                        end
    % error line:             x2=j+l+1;		change to:
                            x2 = j + l;		% The same as above.
                        if x2<1
    % error line              x2=x2-2*l+1;		change to:
                            x2 = -x2 + 2;		% The same as above.
                        end
                        if x2>n
    % error line:             x2=x2-2*l-1;		change to:
                            x2 = 2*n - x2;		% The same as about.
                        end
                        r(i,j,:)=r(i,j,:)+...
                            im(x1,x2,:)*GuassionSmooth(k+radius+1,l+radius+1);
                    end
                end
            end
        end
        
    otherwise

end

% End timing.
t1 = toc;
disp(['Time consumer: ', num2str(t1), ' s']);

% The code here can be shorten, please keep this in mind. -byJackLee
