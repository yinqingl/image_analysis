% helper functinons for automatic stitching of FISH images
%
% Yinqing Li
% yinqingl@csail.mit.edu
% 2016

function [t,tform] = get_sift_affine_t(img,ref)
%use sift to find an affine transformation
%move the first image so that it aligns with the second image
%first image: img
%second image: tmp

tmp = ref;
%get sift descriptor
[d1, l1]=iat_surf(img);
[d2, l2]=iat_surf(tmp);

%match sift features
% [map, matches, imgInd, tmpInd]=iat_match_features_mex(d1,d2,.7);

%limit maximum number of feature points pairs in calculation
n1 = size(d1,1);
n2 = size(d2,1);
max_n = 10000;
i_s = ceil(max_n^2/n2);
i_s = [1:min(i_s,n1-1):n1];

d1_s = d1;
l1_s = l1;

X1 = [];
X2 = [];

for i = [1:length(i_s)-1],
    if length(i_s)-1 > 1,
        display(sprintf('iat match feature, iteration %d/%d',i,length(i_s)-1));
    end
    d1 = d1_s(i_s(i):i_s(i+1),:);
    l1 = l1_s(i_s(i):i_s(i+1),:);
    [~, ~, imgInd, tmpInd]=iat_match_features_mex(d1,d2,.7);
    X1 = [X1;l1(imgInd,1:2)];
    X2 = [X2;l2(tmpInd,1:2)];
end

% [~, ~, imgInd, tmpInd]=iat_match_features_mex(d1,d2,.7);

% X1 = l1(imgInd,1:2);
% X2 = l2(tmpInd,1:2);

%plot correspondences
%{
X1 = [];
X2 = [];
[~, ~, imgInd, tmpInd]=iat_match_features_mex(d1,d2,.999);
X1 = [X1;l1(imgInd,1:2)];
X2 = [X2;l2(tmpInd,1:2)];
iat_plot_correspondences(img, tmp, X1(1:end,:)', X2(1:end,:)');
%}

X1h = iat_homogeneous_coords(X1');
X2h = iat_homogeneous_coords(X2');

%use ransac to remove outlier
% [inliers, ransacWarp]=iat_ransac( X2h, X1h,'affine','tol',.05, 'maxInvalidCount', 100);
w = warning ('on','all');
[inliers, ransacWarp]=iat_ransac( X2h, X1h,'affine','tol',.05, 'maxInvalidCount', 10);
warning(w)

if length(inliers) <= 20,
    error('not enough correspondences points...');
end

%plot correspondences
% iat_plot_correspondences(img,tmp,X1(inliers,:)',X2(inliers,:)');

t = affine2d(inv(ransacWarp)');
tform = maketform('affine',t.T);

%{
%perform affine transformation to register image
[M,N] = size(tmp);
% [wimage, support] = iat_inverse_warping(img, ransacWarp, 'affine', 1:N, 1:M);

%show images side by side
% figure, imshowpair(tmp, wimage, 'montage')
% figure, imshowpair(tmp, wimage)

Rcb = imref2d(size(tmp)+size(img));

t = affine2d(inv(ransacWarp)');
wimage = imwarp(img,t,'outputview',Rcb);

% figure,imshow(wimage)
%}

end