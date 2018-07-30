function imgn = imrotate2(f, angle, match_paras, varargin)
%IMROTATE2 rotates input image, F, by specified angle, ANGLE. The function
%Optimize with vectorization. The original version comes from internet:
%http://www.cnblogs.com/tiandsp/archive/2012/12/03/2800373.html
%
%   R = IMROTATE2(F, ANGLE, MATCH_PARAS)
%   only uses 'bilinear' way to interpolate new gray
%   value of rotated image after resampling the original image.
%   
%   F should be gray image, only, for now. ANGLE should range from -360 degree
%   to 360 degree.
% 
%   MATCH_PARAS is a structure using to pass parameters among matching functions.
%   The parameters in the structure is introduced as below:
% 
%       match_paras.background_mean         Background_mean; it is the mean intensity of
%                                           the background of image, f, which has been
%                                           calculated in the function, GET_TEMPLATE1.
% 
%       match_paras.P                       Number of layers in the image pyramids.
% 
%       match_paras.layer_size_m            A vector stores numbers of rows
%                                           in each pyramid image. The
%                                           first element is the row of
%                                           original image.
%
%       match_paras.layer_size_n            A vector stores numbers of column
%                                           in each pyramid image. The
%                                           first element is the column of
%                                           original image.
% Example: imgn = imrotate2(f, angle, match_paras, varargin)

% Start timing.
tic;

% Check the number of inputs.
error(nargchk(3, 4, nargin));

% Store the input class.

% Deal with specical inputs, like 360 degree.
if angle == 0
    imgn = f;
    toc;
    return
end

if angle == 360 || angle == -360
    imgn = f;
%     figure, imshow(f, [ ]);
    toc;
    return;
end

if angle == 180 || angle == -180
    imgn = f(end:-1:1, end:-1:1);
%     figure, imshow(imgn, [ ]);
    toc;
    return;
end

if angle == 90 || angle == -270
    % Matrix transposition will only make image mirroring.
    imgn = f';
    imgn = imgn(1:end, end:-1:1);
%     figure, imshow(imgn, [ ]);
    toc;
    return;
end

if angle == -90 || angle == 270
    imgn = f';
    imgn = imgn(end:-1:1, 1:end);
%     figure, imshow(imgn, [ ]);
    toc;
    return;
end
    
% Change input calss into 'double', which is hope for MATLAB.
f = double(f);

% Set default values.
% background_mean = 256/2;
% if nargin > 3
%     match_paras = varargin{1};
% end

% Predefine lager varibles.

% Only deal with gray image for now.
[~, ~, imL] = size(f);
if imL ~= 1
    error('The input image should be gray image for now.--imrotate2');
end

% Verify the input angle is within the scope of [-360, 360].
if angle < -360 || angle > 360
    error('The input angle is within the scope of [0, 360].--imrotate2');
end

jiaodu=angle;                       %要旋转的角度，旋转方向为顺时针
img = f;                            %这里v为原图像的高度，u为原图像的宽度
% figure, imshow(img, [ ]);         %这里y为变换后图像的高度，x为变换后图像的宽度
[h w]=size(img);

theta=jiaodu/180*pi;

% A simple way to change the positive rotating directio, here.
% theta = -theta;

% The difference in the expression from the book make switching of 'w' and
% 'h' later. This is the expression in the book:
rot=[cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0;0 0 1];

% Positive rotating direction is anticlockwise:
% rot=[cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0;0 0 1];
% Positive rotating direction is clockwise:
% rot=[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];       

% The positions of ROW and COLUMN have been changed.
pix1=[1 1 1]*rot;               %变换后图像左上点的坐标
pix2=[1 w 1]*rot;               %变换后图像右上点的坐标
pix3=[h 1 1]*rot;               %变换后图像左下点的坐标
pix4=[h w 1]*rot;               %变换后图像右下点的坐标

% Calculate distance and round upwards.
height=round(max([abs(pix1(1)-pix4(1))+0.5 abs(pix2(1)-pix3(1))+0.5]));     %变换后图像的高度
width=round(max([abs(pix1(2)-pix4(2))+0.5 abs(pix2(2)-pix3(2))+0.5]));      %变换后图像的宽度

% Set background color outside image boundary.
% imgn=zeros(height,width);

% Set white background to make wrong edge evident.
% imgn=255*ones(height,width);

% Make vectorizing image.
D = height*width;
imgv = match_paras.background_mean*ones(D, 1);

% delta_y=abs(min([pix1(1) pix2(1) pix3(1) pix4(1)]));            %取得y方向的负轴超出的偏移量
% delta_x=abs(min([pix1(2) pix2(2) pix3(2) pix4(2)]));            %取得x方向的负轴超出的偏移量

rmin = min([pix1(1) pix2(1) pix3(1) pix4(1)]);                    % Range of rotated image.
rmax = max([pix1(1) pix2(1) pix3(1) pix4(1)]);
cmin = min([pix1(2) pix2(2) pix3(2) pix4(2)]);
cmax = max([pix1(2) pix2(2) pix3(2) pix4(2)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%% former FOR loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i = rmin : rmax                           % No sue for(rmax + 1), becouse line 100: 'if pix(i, 1)>=1 && pix(i, 2)>=1 && pix(i, 1) <= h && pix(i, 2) <= w'
%                                               % will let the extra rows or columns pass.
%                                               % There should be a '0' index! Should '0' index be a question.
%     for j = cmin : cmax
% 
% % for i=1-delta_y:height-delta_y              % Why should be like this???
% %     for j=1-delta_x:width-delta_x           % The black edge should due to this. When there is no
%                                             % negative offset, the
%                                             % algorithm here is wrong. Many
%                                             % times, the minimum
%                                             % coordinate of the retated
%                                             % image is 1. But the
%                                             % expression above(delta_y and delta_x) will treat
%                                             % it as '-1'.
%                                             % Confirmed! Your guessing is
%                                             % right. But there still left
%                                             % one single black edge.
%                                             
%         pix=[i j 1]/rot;                                  %用变换后图像的点的坐标去寻找原图像点的坐标，                                         
                                                            %否则有些变换后的图像的像素点无法完全填充
                                                        
%%%%%%%%%%%%%%%%%%%%%%%% Vectorization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[C, R] = meshgrid(cmin:cmax, rmin:rmax);
% To match with matrix 'rot', change the order of R and C.
CR = [R(:), C(:), ones(D, 1)];      
pix = CR/rot;

%%%%%%%%%%%%%%%%%%%%%%%% Padding %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgm=img_extend(img,1);

% Extend original image range threshold a little to deal with rounding error.                                                                       
h2 = h + 0.0001;
w2 = w + 0.0001;

% Use integer as index here. It is better.
for i = 1:D          
% Decrease any calculation in FOR loop, if can.

% %       pix = pix + 0.000001;                     % Offset the rounding error of the inverse operation.
%         if abs(pix(i, 1)-1) < 0.00000000001       % Checking the error every is not good idea. Use padding
%                                                   % will be more efficient.
%             pix(i, 1)=1;
%         end
%         if abs(pix(i, 1)-h) < 0.00000000001 
%             pix(i, 1)=h;
%         end            
%         if abs(pix(i, 2)-1) < 0.00000000001 
%             pix(i, 2)=1;
%         end        
%         if abs(pix(i, 2)-w) < 0.00000000001 
% %           pix(i, 1)=w;                        % Writing error!!! Error typing.
%             pix(i, 2) = w;
%         end        

%       if pix(i, 1)>=1 && pix(i, 2)>=1 && pix(i, 1) <= h && pix(i, 2) <= w

        % Extend range a little to deal with rounding error.
        if pix(i, 1)>=0.9999 && pix(i, 2)>=0.9999 && pix(i, 1) <= h2 && pix(i, 2) <= w2     
                                  
        float_Y=pix(i, 1)-floor(pix(i, 1)); 
        float_X=pix(i, 2)-floor(pix(i, 2));
        
%       if pix(i, 1)>=1 && pix(i, 2)>=1 && pix(i, 1) < h + 1 && pix(i, 2) < w + 1 
%       if round(pix(i, 1))>=1 && round(pix(i, 2))>=1 && pix(i, 1) <= h && pix(i, 2) <= w
%       if pix(i, 1)>0 && pix(i, 2)>0 && pix(i, 1) <= h && pix(i, 2) <= w
                                            % Look like no problem. The rounded decimals hide it.
                                            % '1.00000' in not '1', it is
                                            % '0.99999999999'. So be
                                            % careful, when using decimals
                                            % as counters.

        pix_up_left=[floor(pix(i, 1)) floor(pix(i, 2))];          %四个相邻的点
        pix_up_right=[floor(pix(i, 1)) ceil(pix(i, 2))];
        pix_down_left=[ceil(pix(i, 1)) floor(pix(i, 2))];
        pix_down_right=[ceil(pix(i, 1)) ceil(pix(i, 2))]; 

        value_up_left=(1-float_X)*(1-float_Y);              %计算临近四个点的权重
        value_up_right=float_X*(1-float_Y);
        value_down_left=(1-float_X)*float_Y;
        value_down_right=float_X*float_Y;
        
        % IMGM has been extend to [m+2, n+2]. So add '1' to old index.                                
        imgv(i)=value_up_left*imgm(pix_up_left(1)+1, pix_up_left(2)+1)+ ...          
                value_up_right*imgm(pix_up_right(1)+1,pix_up_right(2)+1)+ ...
                value_down_left*imgm(pix_down_left(1)+1,pix_down_left(2)+1)+ ...
                value_down_right*imgm(pix_down_right(1)+1,pix_down_right(2)+1);
        end       
end

% Inverse vectorization.
imgn = reshape(imgv(:), height, width);
disp(['The original size of rotated image is ', num2str(size(f))]);
disp(['The new size of rotated image is ', num2str(size(imgn))]);
% size(imgn)
% figure,imshow(uint8(imgn), [ ]);
% figure,imshow(imgn, [ ]);

% figure, imshow(im2uint8(mat2gray(imgn)));

disp('Function: imrotate2.');
toc;
% There is another FUNCTION in the M file. So each FUNCTION need a 'END' at
% the end.
end

%%%%%%%%%%%%%%%%%%%%%%%% img_extend %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The extension is similar to 'replicate' option.
function imgm=img_extend(img,r)
[m n]=size(img);

imgm=zeros(m+2*r,n+2*r);

imgm(1+r:m+r,1+r:n+r)=img;
imgm(1:r,1+r:n+r)=repmat(img(1,:), r, 1); 
imgm(1:m+r,n+r+1:n+2*r)=repmat(imgm(1:m+r,n+r), 1, r);
imgm(m+r+1:m+2*r,r:n+2*r)=repmat(imgm(m+r,r:n+2*r), r, 1);
imgm(1:m+2*r,1:r)=repmat(imgm(1:m+2*r,r+1), 1, r);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 90 degree %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r = imrotate2(g, 90);
% ans =
%          800        1024
% Elapsed time is 5.855796 seconds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 90 degree %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imrotate_test3(g, 90);
% ans =
%          800        1024
% Elapsed time is 6.491376 seconds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 90 degree %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% IMROTATE has speed up +/- 90 rotating separately.
% tic;
% r = imrotate(g, -90);
% size(r)
% figure, imshow(r, [ ]);
% toc;
% ans =
%          800        1024
% Elapsed time is 0.294083 seconds.
% IMROTATE use C Language code--'imrotatemex', and may use matrix
% decomposition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 37 degree %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r = imrotate2(g, 37);
% ans =
%         1298        1254
% Elapsed time is 6.238451 seconds.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 37 degree %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imrotate_test3(g, 37);
% ans =
%         1298        1254
% Elapsed time is 8.803251 seconds.        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 37 degree %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic;
% r = imrotate(g, -37);
% figure, imshow(r, [ ]);
% toc;
% Elapsed time is 0.344968 seconds.