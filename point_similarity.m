function s_point = point_similarity(im, im_padded, w, i, j, varargin)
%POINT_SIMILARITY Calculates a point of template matched-degree of specified
%coordinate, (I, J), in the image, IM!!! W is the template. 
%POINT_SIMILARITY transforms from function CALCU_SIMILARITY.
%
%   S_POINT = POINT_SIMILARITY(IM, IM_PADDED, W, SI, SJ). The input image 
%   should be gray. The image must be padded with function PAD_IMAGE.
%   Template, W, should be smaller then input image.
%
%   S = POINT_SIMILARITY(IM, IM_PADDED, W, SI, SJ, METRIC_OPTION). METRIC_OPTION is a string 
%   that can specify a type of similarity metric from below. The default
%   value is enclosed in braces ({}).
% 
%       {'SAD'} Sum of absolute differences.
% 
%       'ZSAD'  Zero-mean sum of absolute differences.
% 
%       'SSD'   Sum of squared differences.
% 
%       'ZSSD'  Zero-mean sum of squared differences.
% 
%       'NCC'   Normalized cross correlation.
% 
%       'ZNCC'  Zero-mean normalized cross correlation.
%       (To see the formulas about metrics, refer to the website mentioned
%       above.)
%
%   S_POINT = POINT_SIMILARITY(IM, IM_PADDED, W,  SI, SJ, METRIC_OPTION, CALCU_METHOD).
%   CALCU_METHOD is a string that can choose from below options to define a
%   type of algorithm calculation stratagem. The default value for
%   CALCU_METHOD is enclosed in braces ({}).
%
%       {'time-save'}   The algorithm will be calculated in a way of
%                       time-save and space-consuming.
%
%       'space-save'    The algorithm will be calculated in a way of
%                       space-save and time-consuming.
%

% Start timing.
% tic;

% Verify the right number of inputs.
error(nargchk(5, 7, nargin));

% Store the format of input image for use later.
classin1 = class(im);
classin2 = class(w);

% Transform class of input image to 'doublt' that MATLAB need for
% calculation.
im = double(im);
w = double(w);

% Set default values.
metric_option = 'SAD';
calcu_method = 'time_save';

% Predefine for large variables and get input values.
if nargin > 6
    calcu_method = varargin{2};
end
if nargin > 5
    metric_option = varargin{1};
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

m_mid = ceil(m/2);
n_mid = ceil(n/2);

padding_h = floor(m/2);
padding_w = floor(n/2);

im_mean = sum(im(:))/NUM;
w_mean = sum(w(:))/num;

sub_im = zeros(m, n);

% sub_mean_im is a array that has the same size of im and whose elements
% all equal the mean of the coorisponding sub-image covered by template.
sub_mean_im = zeros(m, n);
% mean_w is template-size matrix formed by mean of template, w.
mean_w = w_mean*ones(m, n);

%%%%%%%%%%%%%%%%%%%%%% Calculate metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch metric_option
    case 'SAD'
        switch calcu_method
            case 'space_save'
%                 tic;
                sub_im = im_padded(i+padding_h-m_mid+1:i+padding_h-m_mid+m,...
                    j+padding_w-n_mid+1:j+padding_w-n_mid+n);
                s_point = sum(sum(abs(w - sub_im)));
%                 disp('Function: POINT_SIMILARITY; time for calculating SAD, using space-save way.');
%                 toc;
%%%%%%%%%%%%%%%%%%% speed test code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'time_save'
%                 tic;
                m1 = m_mid-1-padding_h;
                m2 = m_mid-m-padding_h;
                n1 = n_mid-1-padding_w;
                n2 = n_mid-n-padding_w;        

                sub_im = im_padded(i-m1:i-m2, j-n1:j-n2);
                s_point = sum(sum(abs(w - sub_im)));
%                 disp('Function: POINT_SIMILARITY; time for calculating SAD, using time-save way.');
%                 toc;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    case 'ZSAD'
        switch calcu_method
            case 'space_save'        
%                 tic;
                sub_im = im_padded(i+padding_h-m_mid+1:i+padding_h-m_mid+m,...
                    j+padding_w-n_mid+1:j+padding_w-n_mid+n);
                sub_mean_im = sum(sub_im(:))/num;
                s_point = sum(sum(abs(w - mean_w - sub_im + sub_mean_im))); 
%                 disp('Function: POINT_SIMILARITY; time for calculating ZSAD, using space-save way.');
%                 toc;        
%%%%%%%%%%%%%%%%%%% speed test code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'time_save'
%                 tic;
                m1 = m_mid-1-padding_h;
                m2 = m_mid-m-padding_h;
                n1 = n_mid-1-padding_w;
                n2 = n_mid-n-padding_w;        
                w3 = w - mean_w;

                sub_im = im_padded(i-m1:i-m2, j-n1:j-n2);
                sub_mean_im = sum(sub_im(:))/num;
                s_point = sum(sum(abs(w3 - sub_im + sub_mean_im))); 
%                 disp('Function: POINT_SIMILARITY; time for calculating ZSAD, using time-save way.');
%                 toc;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'SSD'
        switch calcu_method
            case 'space_save'        
                sub_im = im_padded(i+padding_h-m_mid+1:i+padding_h-m_mid+m,...
                    j+padding_w-n_mid+1:j+padding_w-n_mid+n);
                s_point = sum(sum((w - sub_im).*(w - sub_im)));
                disp('Function: POINT_SIMILARITY; time for calculating SSD, using space-save way.');
%                 toc
%%%%%%%%%%%%%%%%%%% speed test code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'time_save'
                m1 = m_mid-1-padding_h;
                m2 = m_mid-m-padding_h;
                n1 = n_mid-1-padding_w;
                n2 = n_mid-n-padding_w; 
                % Predefine variables outside the FOR loop to accelerate
                % the program. Here this operation may contribute little,
                % however, this is a notice.
                dif_w_im = zeros(m, n); 
                sub_im = im_padded(i-m1:i-m2, j-n1:j-n2);
                dif_w_im = w - sub_im;
                s_point = sum(sum(dif_w_im.*dif_w_im));
%                 disp('Function: POINT_SIMILARITY; time for calculating SSD, using time-save way.');
%                 toc;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'ZSSD'
        switch calcu_method
            case 'space_save'        
                sub_im = im_padded(i+padding_h-m_mid+1:i+padding_h-m_mid+m,...
                    j+padding_w-n_mid+1:j+padding_w-n_mid+n);
                sub_mean_im = sum(sub_im(:))/num;
                s_point = sum(sum((w - mean_w - sub_im + sub_mean_im)...
                    .*(w - mean_w - sub_im + sub_mean_im))); 
                disp('Function: POINT_SIMILARITY; time for calculating ZSSD, using space-save way.');
%                 toc;        
%%%%%%%%%%%%%%%%%%% speed test code %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            case 'time_save'
                m1 = m_mid-1-padding_h;
                m2 = m_mid-m-padding_h;
                n1 = n_mid-1-padding_w;
                n2 = n_mid-n-padding_w; 
                zmean_w = w - mean_w;
                
                sub_im = im_padded(i-m1:i-m2, j-n1:j-n2);
                sub_mean_im = sum(sub_im(:))/num;
                % Error code!!! Be careful to accumulatiom in FOR
                % loops, like below: 'zmean_zssd = zmean_zssd...'.
                % zmean_zssd should not be accumulation here. The
                % way of calculation zmean_zssd will accumulate its 
                % result of last calculation. Error!!!
                % 'zmean_zssd = zmean_zssd - sub_im + sub_mean_im;'
                zmean_zssd = zmean_w - sub_im + sub_mean_im;
                s_point = sum(sum(zmean_zssd.*zmean_zssd)); 
%                 disp('Function: POINT_SIMILARITY; time for calculating ZSSD, using time-save way.');
%                 toc;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'NCC'
        switch calcu_method
            case 'space_save'        
                sub_im = im_padded(i+padding_h-m_mid+1:i+padding_h-m_mid+m,...
                    j+padding_w-n_mid+1:j+padding_w-n_mid+n);
                s_point = sum(sum(w.*sub_im))/(sqrt(sum(sum(w.*w))*sum(sum(sub_im.*sub_im))));
                disp('Function: POINT_SIMILARITY; time for calculating NCC, using space-save way.');
%                 toc;        
%%%%%%%%%%%%%%%%%%% speed test code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'time_save'
                m1 = m_mid-1-padding_h;
                m2 = m_mid-m-padding_h;
                n1 = n_mid-1-padding_w;
                n2 = n_mid-n-padding_w; 
                ww_ncc = sum(sum(w.*w));
                
                sub_im = im_padded(i-m1:i-m2, j-n1:j-n2);
                s_point = sum(sum(w.*sub_im))/(sqrt(ww_ncc*sum(sum(sub_im.*sub_im))));
%                 disp('Function: POINT_SIMILARITY; time for calculating NCC, using time-save way.');
%                 toc;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'ZNCC'
        switch calcu_method
            case 'space_save'        
                sub_im = im_padded(i+padding_h-m_mid+1:i+padding_h-m_mid+m,...
                    j+padding_w-n_mid+1:j+padding_w-n_mid+n);
                sub_mean_im = sum(sub_im(:))/num;
                s_point = sum(sum((w-mean_w).*(sub_im-sub_mean_im)))/...
                    sqrt(sum(sum((w-mean_w).*(w-mean_w)))*sum(sum((sub_im-sub_mean_im).*(sub_im-sub_mean_im))));
                disp('Function: POINT_SIMILARITY; time for calculating ZNCC, using space-save way.');
%                 toc;  
%%%%%%%%%%%%%%%%%%% speed test code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'time_save'
%                 tic;
                m1 = m_mid-1-padding_h;
                m2 = m_mid-m-padding_h;
                n1 = n_mid-1-padding_w;
                n2 = n_mid-n-padding_w;
                norm_w = w-mean_w;
                norm_ww = norm_w.*norm_w;

                sub_im = im_padded(i-m1:i-m2, j-n1:j-n2);
                sub_mean_im = sum(sub_im(:))/num;
                norm_sub_im = sub_im-sub_mean_im;
                s_point = sum(sum(norm_w.*norm_sub_im))/...
                sqrt(sum(sum(norm_ww))*sum(sum(norm_sub_im.*norm_sub_im)));
%                 disp('Function: POINT_SIMILARITY; time for calculating ZNCC, using time-save way.');
%                 toc;        
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    otherwise
        error('Unknown similarity metric has been input.(POINT_SIMILARITY)');
end

% Output.

% disp('Function: POINT_SIMILARITY; total elapsed time.');
% % End timing.
% toc;




        