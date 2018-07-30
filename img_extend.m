function imgm=img_extend(img,r)
%IMG_EXTEND Pads the image, IMG, by copy edge pixel R times.
%IMGM = IMG_EXTEND(IMG, R)

[m n]=size(img);

imgm=zeros(m+2*r,n+2*r);

imgm(1+r:m+r,1+r:n+r)=img;
imgm(1:r,1+r:n+r)=repmat(img(1,:), r, 1); 
imgm(1:m+r,n+r:n+2*r)=repmat(imgm(1:m+r,n+r), 1, r);
imgm(m+r+1:m+2*r,r:n+2*r)=repmat(imgm(m+r,r:n+2*r), r, 1);
imgm(1:m+2*r,1:r)=repmat(imgm(1:m+2*r,r+1), 1, r);

end