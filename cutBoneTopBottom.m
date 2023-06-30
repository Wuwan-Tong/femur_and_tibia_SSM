function [cutBoneSet] = cutBoneTopBottom(LocMatrixSet,z_max_perc,z_min_perc)
% cut off the top part and the bottom part of the bone
% input: 
%    LocMatrixSet: 1-by-mun_of_img cell array, each cell is the location matrix of one point cloud image
%    z_max_perc,z_min_perc: double, in percent, e.g. if z_max=0.6, z_min=0.1, then cut 40% of the bone at top (z-axis) and 10% at bottom
% output: 
%    cutBoneSet: 1-by-mun_of_img cell array, each cell stores the location matrix of one point cloud image after cutting

numImage=length(LocMatrixSet); %the number of input images
% create an empty downSampBone
cutBoneSet=cell(1,numImage);

for i=1:numImage
    tempBone=cell2mat(LocMatrixSet(i));
    % find z_max and z_min according to z_max_percent and z_min_percent for each bone
    z_max=(max(tempBone(:,3))-min(tempBone(:,3)))*z_max_perc+min(tempBone(:,3));
    z_min=(max(tempBone(:,3))-min(tempBone(:,3)))*z_min_perc+min(tempBone(:,3));
    % cut 
    tempBone(find(tempBone(:,3)>z_max|tempBone(:,3)<z_min),:)=[];
    % store in cutBoneSet
    [rows,cols]=size(tempBone);
    cutBoneSet(1,i)=mat2cell(tempBone,[rows],[cols]);
end
end