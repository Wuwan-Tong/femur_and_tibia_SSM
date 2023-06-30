function [cutBone] = cut_top(LocMatrixSet,cutOffPerc,toScale)
% cut the transformed femur point cloud after TMM, cut the z>z_max part

% Input: 
%    LocMatrixSet: 1-by-K cell array, location matrix of point clouds
%    z_max, cut z>z_max part of the femur

% Output:
%    cutBone, 1-by-k cell array of n-by-3 point cloud location, with n<N, the
%    number of rest points after cutting

K=length(LocMatrixSet);
% create an empty cutBone
cutBone=cell(1,K);
figure
if toScale
    %find z max
    tempBone=cell2mat(LocMatrixSet(1));
    z_max=(max(tempBone(:,3))-min(tempBone(:,3)))*(1-cutOffPerc)+min(tempBone(:,3));
    % cut    
    for i=1:K
        tempBone=cell2mat(LocMatrixSet(i));
        % scatter3(tempBone(:,1),tempBone(:,2),tempBone(:,3),'.');
        tempBone(find(tempBone(:,3)>z_max),:)=[];
        [rows,cols]=size(tempBone);
        cutBone(1,i)=mat2cell(tempBone,[rows],[cols]);
        scatter3(tempBone(:,1),tempBone(:,2),tempBone(:,3),'.');
        hold on
    end
    title('bones after prealignment and cut')
else
    tempBone=cell2mat(LocMatrixSet(1));
    z_min_first_bone=min(tempBone(:,3));
    z_length=(max(tempBone(:,3))-z_min_first_bone)*(1-cutOffPerc); % the length of bone left
    for i=1:K
        tempBone=cell2mat(LocMatrixSet(i));
        %find z max
        z_min=min(tempBone(:,3));
        z_max=z_length+z_min;
    %     scatter3(tempBone(:,1),tempBone(:,2),tempBone(:,3),'.');
        tempBone(find(tempBone(:,3)>z_max),:)=[];
        [rows,cols]=size(tempBone);
        cutBone(1,i)=mat2cell(tempBone,[rows],[cols]);
        % show bones move to same z_min
        tempBone_print=tempBone;
        tempBone_print(:,3)=tempBone(:,3)-(z_min-z_min_first_bone);
        scatter3(tempBone_print(:,1),tempBone_print(:,2),tempBone_print(:,3),'.');
        hold on
    end
    title('bones after prealignment and cut, and aligned again along z')
end
daspect([1 1 1]);

