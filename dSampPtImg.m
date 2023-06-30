function [dSampBone] = dSampPtImg(LocMatrixSet,dsamDist)
% downsampling the point cloud images
% input: 
%    LocMatrixSet: 1-by-mun_of_img cell array, each cell is the location matrix of one point cloud image
%    dsamDist: double, downsampling distance
% output: 
%    dSampBone: 1-by-mun_of_img cell array, each cell is the location matrix of downsampled point cloud image
numImage=length(LocMatrixSet); %the number of input images
% create an empty downSampBone
dSampBone=cell(1,numImage);
for i=1:numImage
    % switch location matrix to point cloud
    pt=pointCloud(cell2mat(LocMatrixSet(i)));
    % downsampling
    dsampPt=pcdownsample(pt,'gridAverage',dsamDist);
    % store the location of downsampled image
    [rows,cols]=size(dsampPt.Location);
    dSampBone(i)=mat2cell(dsampPt.Location,[rows],[cols]);
end
end