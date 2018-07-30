function [s, target_num, I2, J2, peak_values2, polarity, match_info] = calcu_similarity2(im, w, multi_mode, threshold_pecentage, match_paras, varargin)
%CALCU_SIMILARITY Calculates similarities of template, W, at different
%location in the target image, IM, using specified metric of similarity. To
%get more information about different metrics of similarity, refer to
%internet website in Chrome:
%https://siddhantahuja.wordpress.com/tag/sum-of-squared-differences/.
%
%   S = CALCU_SIMILARITY(IM, W). The input image should be gray.
%   Template, W, should be smaller then input image.
%
%   S = CALCU_SIMILARITY(IM, W, METRIC_OPTION). METRIC_OPTION is a string 
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
%   S = CALCU_SIMILARITY(IM, W, METRIC_OPTION, CALCU_METHOD).
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
%   S = CALCU_SIMILARITY(IM, W, METRIC_OPTION, CALCU_METHOD, PADDING_MATHOD, B).
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
error(nargchk(5, 9, nargin));

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
padding_method = 'replicate';
b = 0;

% Predefine for large variables and get input values.
if nargin > 8
    b = varargin{4};
end
if nargin > 7
    padding_method = varargin{3};
end
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

im_padded = 255*ones(M + 2*padding_h, N + 2*padding_w);
im_padded(padding_h+1:padding_h+M, padding_w+1:padding_w+N) = im;

im_mean = sum(im(:))/NUM;
w_mean = sum(w(:))/num;

s = zeros(M, N);
sub_im = zeros(m, n);

% sub_mean_im is a array that has the same size of im and whose elements
% all equal the mean of the coorisponding sub-image covered by template.
sub_mean_im = zeros(m, n);
% mean_w is template-size matrix formed by mean of template, w.
mean_w = w_mean*ones(m, n);

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
        error('Sorry! The symmetric padding function is not finished.');
        
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
        error('Unknown padding method has been inputed.(calcu_similarity)');
end

%%%%%%%%%%%%%%%%%%%%%% Calculate metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch metric_option
    case 'SAD'
        switch calcu_method
            case 'space_save'
                tic;
                for i = 1:M
                    for j = 1:N
                        sub_im = im_padded(i+padding_h-m_mid+1:i+padding_h-m_mid+m,...
                            j+padding_w-n_mid+1:j+padding_w-n_mid+n);
                        s(i, j) = sum(sum(abs(w - sub_im)));
                    end
                end
                disp('Function: calcu_similarity; time for calculating SAD, using space-save way.');
                toc;
%%%%%%%%%%%%%%%%%%% speed test code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'time_save'
                tic;
                m1 = m_mid-1-padding_h;
                m2 = m_mid-m-padding_h;
                n1 = n_mid-1-padding_w;
                n2 = n_mid-n-padding_w;        

                for i = 1:M
                    for j = 1:N
                        sub_im = im_padded(i-m1:i-m2, j-n1:j-n2);
                        s(i, j) = sum(sum(abs(w - sub_im)));
                    end
                end
                disp('Function: calcu_similarity; time for calculating SAD, using time-save way.');
                toc;
        end
        [I, J] = find(s == min(s(:)));
        polarity = 'minimum';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    case 'ZSAD'
        switch calcu_method
            case 'space_save'        
                tic;
                for i = 1:M
                    for j = 1:N
                        sub_im = im_padded(i+padding_h-m_mid+1:i+padding_h-m_mid+m,...
                            j+padding_w-n_mid+1:j+padding_w-n_mid+n);
                        sub_mean_im = sum(sub_im(:))/num;
                        s(i, j) = sum(sum(abs(w - mean_w - sub_im + sub_mean_im)));        
                    end
                end
                disp('Function: calcu_similarity; time for calculating ZSAD, using space-save way.');
                toc;        
%%%%%%%%%%%%%%%%%%% speed test code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'time_save'
                tic;
                m1 = m_mid-1-padding_h;
                m2 = m_mid-m-padding_h;
                n1 = n_mid-1-padding_w;
                n2 = n_mid-n-padding_w;        
                w3 = w - mean_w;

                for i = 1:M
                    for j = 1:N
                        sub_im = im_padded(i-m1:i-m2, j-n1:j-n2);
                        sub_mean_im = sum(sub_im(:))/num;
                        s(i, j) = sum(sum(abs(w3 - sub_im + sub_mean_im)));        
                    end
                end
                disp('Function: calcu_similarity; time for calculating ZSAD, using time-save way.');
                toc;
        end
        [I, J] = find(s == min(s(:)));
        polarity = 'minimum';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'SSD'
        switch calcu_method
            case 'space_save'        
                for i = 1:M
                    for j = 1:N
                        sub_im = im_padded(i+padding_h-m_mid+1:i+padding_h-m_mid+m,...
                            j+padding_w-n_mid+1:j+padding_w-n_mid+n);
                        s(i, j) = sum(sum((w - sub_im).*(w - sub_im)));        
                    end
                end
                disp('Function: calcu_similarity; time for calculating SSD, using space-save way.');
                toc; 
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
                for i = 1:M
                    for j = 1:N
                        sub_im = im_padded(i-m1:i-m2, j-n1:j-n2);
                        dif_w_im = w - sub_im;
                        s(i, j) = sum(sum(dif_w_im.*dif_w_im));        
                    end
                end
                disp('Function: calcu_similarity; time for calculating SSD, using time-save way.');
                toc;
        end
        [I, J] = find(s == min(s(:)));
        polarity = 'minimum';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'ZSSD'
        switch calcu_method
            case 'space_save'        
                for i = 1:M
                    for j = 1:N
                        sub_im = im_padded(i+padding_h-m_mid+1:i+padding_h-m_mid+m,...
                            j+padding_w-n_mid+1:j+padding_w-n_mid+n);
                        sub_mean_im = sum(sub_im(:))/num;
                        s(i, j) = sum(sum((w - mean_w - sub_im + sub_mean_im)...
                            .*(w - mean_w - sub_im + sub_mean_im)));             
                    end
                end
                disp('Function: calcu_similarity; time for calculating ZSSD, using space-save way.');
                toc;        
%%%%%%%%%%%%%%%%%%% speed test code %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            case 'time_save'
                m1 = m_mid-1-padding_h;
                m2 = m_mid-m-padding_h;
                n1 = n_mid-1-padding_w;
                n2 = n_mid-n-padding_w; 
                zmean_w = w - mean_w;
                for i = 1:M
                    for j = 1:N
                        sub_im = im_padded(i-m1:i-m2, j-n1:j-n2);
                        sub_mean_im = sum(sub_im(:))/num;
                        % Error code!!! Be careful to accumulatiom in FOR
                        % loops, like below: 'zmean_zssd = zmean_zssd...'.
                        % zmean_zssd should not be accumulation here. The
                        % way of calculation zmean_zssd will accumulate its 
                        % result of last calculation. Error!!!
                        % 'zmean_zssd = zmean_zssd - sub_im + sub_mean_im;'
                        zmean_zssd = zmean_w - sub_im + sub_mean_im;
                        s(i, j) = sum(sum(zmean_zssd.*zmean_zssd));       
                    end
                end
                disp('Function: calcu_similarity; time for calculating ZSSD, using time-save way.');
                toc;
        end
        [I, J] = find(s == min(s(:)));
        polarity = 'minimum';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'NCC'
        switch calcu_method
            case 'space_save'        
                for i = 1:M
                    for j = 1:N
                        sub_im = im_padded(i+padding_h-m_mid+1:i+padding_h-m_mid+m,...
                            j+padding_w-n_mid+1:j+padding_w-n_mid+n);
                        s(i, j) = sum(sum(w.*sub_im))/(sqrt(sum(sum(w.*w))*sum(sum(sub_im.*sub_im))));
                    end
                end
                disp('Function: calcu_similarity; time for calculating NCC, using space-save way.');
                toc;        
%%%%%%%%%%%%%%%%%%% speed test code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'time_save'
                m1 = m_mid-1-padding_h;
                m2 = m_mid-m-padding_h;
                n1 = n_mid-1-padding_w;
                n2 = n_mid-n-padding_w; 
                ww_ncc = sum(sum(w.*w));
                for i = 1:M
                    for j = 1:N
                        sub_im = im_padded(i-m1:i-m2, j-n1:j-n2);
                        s(i, j) = sum(sum(w.*sub_im))/(sqrt(ww_ncc*sum(sum(sub_im.*sub_im))));
                    end
                end 
                disp('Function: calcu_similarity; time for calculating NCC, using time-save way.');
                toc;
        end
        [I, J] = find(s == man(s(:)));
        polarity = 'maximum';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    case 'ZNCC'
        switch calcu_method
            case 'space_save'        
                for i = 1:M
                    for j = 1:N
                        sub_im = im_padded(i+padding_h-m_mid+1:i+padding_h-m_mid+m,...
                            j+padding_w-n_mid+1:j+padding_w-n_mid+n);
                        sub_mean_im = sum(sub_im(:))/num;
                        s(i, j) = sum(sum((w-mean_w).*(sub_im-sub_mean_im)))/...
                            sqrt(sum(sum((w-mean_w).*(w-mean_w)))*sum(sum((sub_im-sub_mean_im).*(sub_im-sub_mean_im))));
                    end
                end
                disp('Function: calcu_similarity; time for calculating ZNCC, using space-save way.');
                toc;  
%%%%%%%%%%%%%%%%%%% speed test code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            case 'time_save'
                tic;
                m1 = m_mid-1-padding_h;
                m2 = m_mid-m-padding_h;
                n1 = n_mid-1-padding_w;
                n2 = n_mid-n-padding_w;
                norm_w = w-mean_w;
                norm_ww = norm_w.*norm_w;

                for i = 1:M
                    for j = 1:N
                        sub_im = im_padded(i-m1:i-m2, j-n1:j-n2);
                        sub_mean_im = sum(sub_im(:))/num;
                        norm_sub_im = sub_im-sub_mean_im;
                        s(i, j) = sum(sum(norm_w.*norm_sub_im))/...
                        sqrt(sum(sum(norm_ww))*sum(sum(norm_sub_im.*norm_sub_im)));
                    end
                end
                disp('Function: calcu_similarity; time for calculating ZNCC, using time-save way.');
                toc;        
        end
        [I, J, peak_values] = find(s == man(s(:)));
        polarity = 'maximum';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    otherwise
        error('Unknown similarity metric has been input.(calcu_similarity)');
end

% Output.
match_info.polarity = polarity;

for sk = 1:length(I)
    peak_values(sk) = s(I(sk), J(sk));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch match_paras.P
    case 1
        threshold_offset = 150000;
    case 2
        threshold_offset = 50000;
    case 3
        threshold_offset = 10000;
    case 4
        threshold_offset = 5000;
    case 5
        threshold_offset = 1000;
    case 6
        threshold_offset = 500;
    otherwise
        error('Sorry, other offsets of higher pyramid laryer is not ready--calcu_similarity2.')
end

switch polarity
    case 'minimum'
        sign = 1;
        s_threshold = (threshold_offset * (1 - threshold_pecentage)) + peak_values(1);
    case 'maximum'
        sign = -1;
        s_threshold = peak_values(1) - (threshold_offset * (1 - threshold_pecentage));
    otherwise
        error('Unknown polarity of similarity.--top_match');
end

switch multi_mode
    case 'off'
        target_num = 1;
        I2 = I;
        J2 = J;
        peak_values2 = peak_values(1);
        
    case 'on'
        target_num = 1;
        I2 = {I(1)};
        J2 = {J(1)};
        peak_values2 = {peak_values(1)};
        
        
        
        for i_s = 1:NUM
            if (s(i_s)*sign) < (s_threshold*sign)
                target_flag = 1;
                for i_target = 1:target_num
                    % If the distance between new matching point and
                    % existed matching target is larger than width of template
                    % witheout rotating, we think the new matching point is a
                    % new target.
                    if sqrt((I2{i_target}-rem(i_s, M))^2 + (J2{i_target}-ceil(i_s/M))^2) < min(size(w))
                        % If distance between new matching point and existed
                        % matching target smaller than width of template,
                        % the two matching points  belongs to the same
                        % target. That comparing their matching similarity,
                        % choose the better.
                        target_flag = 0;
                        % Comparing two matching similarity. Choose better.
                        if (s(i_s)*sign) < (peak_values2{i_target}*sign)
                            I2{i_target} = rem(i_s, M);
                            J2{i_target} = ceil(i_s/M);
                            peak_values2{i_target} = s(i_s);
                            
                        end
                        % Within range of width of template, there can be 
                        % only one target whose distance between new
                        % matching point is smaller than width of
                        % template. End the loop.
                        break
                    end
                end
                % If the new matching point belongs to no exist points, so
                % add it into the target array.
                if target_flag == 1
                    target_num = target_num + 1;
                    I2{target_num} = rem(i_s, M);
                    J2{target_num} = ceil(i_s/M);
                    peak_values2{target_num} = s(i_s);
                end
            end
        end
    otherwise
        error('Unknown option for MULTI_MODE--top_match.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure, imshow(s, [ ]);
% s_t = s';
% figure, surf(s_t(1:round(end/50):end, 1:round(end/50):end));

disp(['There is/are ', num2str(target_num), ' point/points matched.']); 
disp('The matched point(s) and the corresponding peak-metric value(s) is: ');
for ii = 1:target_num
    disp(['(', num2str(I2{ii}), ', ', num2str(J2{ii}), ') ']);
    disp(num2str(peak_values2{ii}));
end

disp('Function: calcu_similarity; total elapsed time.');
% End timing.
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert peak_values2 from cell to array, so can pass it to another cell,
% which can easy index.
peak_values2 = cell2mat(peak_values2);
I2 = cell2mat(I2);
J2 = cell2mat(J2);




        