function [] = plotOneImg(img,title_plot)
% plot one 3D point cloud image using scatter3()
% input: 
%    img: location matrix of the 3D point cloud image
%    title_plot: title of the plot
figure
scatter3(img(:,1),img(:,2),img(:,3),'.');
daspect([1 1 1])
title(title_plot)
end