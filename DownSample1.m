function r=DownSample1(im)
%DOWNSAMPLE1 Performs down-sample of image, IM. Code comes form internet:
%"http://wenku.baidu.com/link?url=RZFBcs_pi
%7apLuafXuBLqMm4VvgL1i4hK7fNWo7DsYIC3flX4b1hdC8DOdPlS894XdeVKzRJgZlE0-Ui
%Gazc7Jx_UKdBI8qhZapWa7cjFwO".
%
% R = DOWNSAMPLE1(IM) downsamples the input image, IM, by nearest neighbor
% method.

	tic;
    [m,n,z]=size(im);
    r=zeros(m/2,n/2,z);
    for i=1:m/2
        for j=1:n/2
            r(i,j,:)=im(2*i,2*j,:);
        end
    end
	toc;
