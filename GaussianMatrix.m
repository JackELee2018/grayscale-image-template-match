function r=GaussianMatrix(delta,radius)
%GAUSSIANMATRIX Creates a gaussian matrix for ImagePyramid algorithm.
%Copyright belongs to Internet "http://wenku.baidu.com/link?url=RZFBcs_pi
%7apLuafXuBLqMm4VvgL1i4hK7fNWo7DsYIC3flX4b1hdC8DOdPlS894XdeVKzRJgZlE0-Ui
%Gazc7Jx_UKdBI8qhZapWa7cjFwO".
%   
%   R = GAUSSIANMATRIX(DELTA, RADIUS) Makes a Approximation Filter for
%   Image Pyramid algorithm. Filtering effect of R is stranger than
%   that of filters that makes by function FSPECIAL.
%   
%   R is output filter matrix.
%   
%   DELTA is sig of Gaussian Expression. DELTA is standard deviation of
%   Gaussian function.
%   
%   RADIUS is radius of the filter matrix.--byJackLee 2016.1.12
 
% Verify the correct number of inputs.--byJackLee 2016.1.12
error(nargchk(2, 2, nargin));

    radius=ceil(radius);
    n=2*radius+1;
    r=zeros(n,n);
    %tempMatrix=zeros(radius,radius);
    for i=-radius:radius
        for j=-radius:radius
            r(i+radius+1,j+radius+1)=exp(-(i^2+j^2)/2*delta^2);
        end
    end;
    r=round(100*r);   
    r=r/sum(sum(r));    % Normalize r to [0, 1].
