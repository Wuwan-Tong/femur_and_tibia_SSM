function [] = plotThreeImg(img1,img2,img3,title_plot)
% plot three 3D point cloud images using scatter3()
% input: 
%    img1,img2,img3: location matrix of the 3D point cloud images
%    title_plot: title of the plot

figure
scatter3(img1(:,1),img1(:,2),img1(:,3),'.');
hold on
scatter3(img2(:,1),img2(:,2),img2(:,3),'.');
scatter3(img3(:,1),img3(:,2),img3(:,3),'.');
daspect([1 1 1])
title(title_plot)
end