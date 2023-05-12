function [] = plotTwoImg(img1,img2,title_plot)
% plot two 3D point cloud images using scatter3()
% input: 
%    img1,img2: location matrix of the 3D point cloud images
%    title_plot: title of the plot

figure
scatter3(img1(:,1),img1(:,2),img1(:,3),'.');
hold on
scatter3(img2(:,1),img2(:,2),img2(:,3),'.');
daspect([1 1 1])
title(title_plot)
end
