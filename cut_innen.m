%% segment innen and outer layer, remove the outliers, then downsample the clean innen and outer layer respectively
% Written by Wuwan Tong 12/05/2023, rwth Aachen

%% set params
% input the dir name where the input images are stored
dirNameTwolayer='aussen_tibia';
dirNameIn='innen_tibia';
% set downsampling distance
dsamDist_femur_out=3.4;
dsamDist_femur_in=2.4;
dsamDist_tibia_out=2.5;
dsamDist_tibia_in=1.7;
% set the dir name where to store the clean downsampled innen and outer
% layer, images are stored in 'dir_path_out/r' and 'dir_path_out/l'
dir_path_out='dSampOut_tibia';
dir_path_in='dSampIn_tibia';
% cut off the top and bottom part, for now we do cutting after prealignment instead, so the 
% code for cutting in this part is commented out, you can reuse this part
% later if it's required

% z_max_perc=0.5; % range:[0,1], the region z>z_max will be cut off, 0.6 means top 40% of the bone will be cut off
% z_min_perc=0.1; % range:[0,1], the region z<z_min will be cut off, 0.1 means bottom 10% of the bone will be cut off
%% content of different cell arrays:
% LocMatrixTwolayer, LocMatrixIn: 
%     cell(1,numImage), store the location of original point clouds
% imgNameSet:
%     cell(1,numImage), store the name of original point clouds, e.g.: 'patient1_femur_r'
% idxOutliers: 
%     cell(1,numImage), store the indices of the outliers(for outer layer) in LocMatrixTwolayer
% LocMatrixOut,LocMatrixIn: 
%     cell(1,numImage), store the location of outer or innen layers, with outliers
% LocMatrixOutClean, LocMatrixInClean : 
%     cell(1,numImage), store the location of outer or innen layers, after trying to remove some outliers
% cutBoneOut,cutBoneIn: 
%     cell(1,numImage), store the location of outer or innen layers, after cutting the top and bottom part
% dSampBoneOut, dSampBoneIn: 
%     cell(1,numImage), store the downsampled cutBoneOut or cutBoneIn

%% Brief intorductions of functions
% function [cutBoneSet] = cutBoneTopBottom(LocMatrixSet,z_max_perc,z_min_perc)
% cut off the top part and the bottom part of the bone
% input: 
%    LocMatrixSet: 1-by-mun_of_img cell array, each cell is the location matrix of one point cloud image
%    z_max_perc,z_min_perc: double, in percent, e.g. if z_max=0.6, z_min=0.1, then cut 40% of the bone at top (z-axis) and 10% at bottom
% output: 
%    cutBoneSet: 1-by-mun_of_img cell array, each cell stores the location matrix of one point cloud image after cutting

% function [dSampBone] = dSampPtImg(LocMatrixSet,dsamDist)
% downsampling the point cloud images
% input: 
%    LocMatrixSet: 1-by-mun_of_img cell array, each cell is the location matrix of one point cloud image
%    dsamDist: double, downsampling distance
% output: 
%    dSampBone: 1-by-mun_of_img cell array, each cell is the location matrix of downsampled point cloud image

% function [] = saveLoc2ply(LocMatrixSet,imgNameSet,dir_path,isInnen)
% downsampling the point cloud images
% input: 
%    LocMatrixSet: 1-by-mun_of_img cell array, each cell is the location matrix of one point cloud image
%    imgNameSet: 1-by-mun_of_img cell array, each cell is the corresponding name of point cloud image
%    dir_path: string, path of dir to save ply images

% function [] = plotOneImg(img,title_plot)
% plot one 3D point cloud image using scatter3()
% input: 
%    img: location matrix of the 3D point cloud image
%    title_plot: title of the plot

% function [] = plotTwoImg(img1,img2,title_plot)
% plot two 3D point cloud images using scatter3()
% input: 
%    img1,img2: location matrix of the 3D point cloud images
%    title_plot: title of the plot

% function [] = plotThreeImg(img1,img2,img3,title_plot)
% plot three 3D point cloud images using scatter3()
% input: 
%    img1,img2,img3: location matrix of the 3D point cloud images
%    title_plot: title of the plot
%% remove the points from 2-layer image if these points also belongs to innen layer image, and get a clean innen layer
% find all images in and both layers
ImgListTwolayer=dir(dirNameTwolayer);
ImgListIn=dir(dirNameIn);
numImage=length(ImgListTwolayer)-2;
% empty cell array to store the Pts locations and names of every images
LocMatrixTwolayer=cell(1,numImage);
LocMatrixIn=cell(1,numImage);
imgNameSet=cell(1,numImage);
% empty cell array to store clean innen layer
LocMatrixInClean=LocMatrixIn;
% create a LocMatrixOut for further calculation by deleting the innen points
LocMatrixOut=LocMatrixTwolayer; 
i=0;% temp variable to count the position to store the Pts location
idxOutliers=cell(1,numImage);% store the index of outliers for each image

for n=1:numImage+2
    i=i+1;
    ImgName=ImgListTwolayer(n).name;
    if ImgName=='.'
        i=i-1;
        continue;
    end
    ImgPatTwolayer=fullfile(dirNameTwolayer,ImgName);
    ImgPathIn=strrep(ImgPatTwolayer,'aussen','innen');
    imgNameSet(i)=cellstr(strrep(ImgName,'_aussen',''));
    % load images to point cloud and store the location in a matrix
    ptCloudTwolayer= pcread(ImgPatTwolayer);
    ptCloudIn= pcread(ImgPathIn);
    LocMatrixTwolayertemp=ptCloudTwolayer.Location; % temp matrix 
    LocMatrixIntemp=ptCloudIn.Location;

% % correct the direction of tibia, 
%     if i~=1&&i~=2&&i~=4&&i~=16&&i~=17&&i~=31 % please recheck the image numbers
%         LocMatrixTwolayertemp(:,3)=-LocMatrixTwolayertemp(:,3);
%         LocMatrixIntemp(:,3)=-LocMatrixIntemp(:,3);
%     end

% %   test start ----------------------------------------------
% display all of the bone in dir, one bone on one figure, Please take
% care, a lot of figures will come to you!!!!!!!!!!!!! and your matlab
% may be down...(@_@)
%     figure
%     scatter3(LocMatrixTwolayertemp(:,1),LocMatrixTwolayertemp(:,2),LocMatrixTwolayertemp(:,3),'.');
%     title(ImgName)
%     daspect([1 1 1])
%     % disp(i) % this number corressponds to the legend at line 157, so
%     you can find which bone in the folder is not as you expected, with
%     printed name by the next line code
%     % disp(ImgName)
% %   test end ---------------------------------------------
    [rows,cols]=size(LocMatrixTwolayertemp);
    LocMatrixTwolayer(i)=mat2cell(LocMatrixTwolayertemp,[rows],[cols]);
    [rows,cols]=size(LocMatrixIntemp); % temp matrix   
    LocMatrixIn(i)=mat2cell(LocMatrixIntemp,[rows],[cols]);
    
    % find the index of points appear in both Innen and Two layers images
    [~,indB] = ismember(cell2mat(LocMatrixIn(i)),cell2mat(LocMatrixTwolayer(i)),'rows'); 
    % get a clean innen layer
    LocMatrixInCleantemp=LocMatrixIntemp(indB~=0,:);
    [rows,cols]=size(LocMatrixInCleantemp);
    LocMatrixInClean(i)=mat2cell(LocMatrixInCleantemp,[rows],[cols]);
    
    countZeros_indB=sum(indB(:)==0);% the number of points that in Innen but not in Both
    temparray=find(indB(:)==0);
    % % test start----------------------------------
    % to show which bone has too many outliers
    % if length(temparray)>14000
    %     disp(length(temparray));
    %     disp(ImgName);
    % end
    % % Test end------------------------------------
    idxOutliers(i)=mat2cell(temparray,length(temparray),1);% find the index of points that in Innen but not in Both    
    indB(cell2mat(idxOutliers(i)))=[];

    % delete the points that appear in both Innen and Two layers images
    LocMatrixOuttemp=cell2mat(LocMatrixTwolayer(i));
    LocMatrixOuttemp(indB,:)=[];
    [rows,cols]=size(LocMatrixOuttemp);
    LocMatrixOut(i)=mat2cell(LocMatrixOuttemp,[rows],[cols]);  
end
% %   test start ----------------------------------------------
% daspect([1 1 1])
% legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35')
% %   test end ---------------------------------------------
%% plot the first image in dir, with innen layer(.), outer layer(.), and outliers(*) of innen layer
temp_img=cell2mat(LocMatrixIn(1));
title_plot='original bone';
plotThreeImg(cell2mat(LocMatrixOut(1)),cell2mat(LocMatrixIn(1)),temp_img(cell2mat(idxOutliers(1)),:),title_plot)
legend('out layer with outliers','clean innen layer','outliers of innen layer')

%% remove the outliers for out layers
tic
LocMatrixOutClean=LocMatrixOut;
for i=1:numImage
    LocMatrixOuttemp=cell2mat(LocMatrixOut(i));
    LocMatrixIntemp=cell2mat(LocMatrixIn(i));
    
    max_dist=0.8; % the size of search space for findNearestNeighbors()
    k=1; % for each outliers of innen layer, find the nearest k point(s) from out layer, and delete them 
    tempidx=cell2mat(idxOutliers(i));
    testLength=zeros(1,length(tempidx));
    idxs=zeros(1,length(tempidx)*k);
    for idx=1:length(tempidx)
        point=LocMatrixIntemp(tempidx(idx),:);
        probPoints=find((LocMatrixOuttemp(:,3)<point(3)+max_dist)&LocMatrixOuttemp(:,2)<point(2)+max_dist&LocMatrixOuttemp(:,1)<point(1)+max_dist&LocMatrixOuttemp(:,3)>point(3)-max_dist&LocMatrixOuttemp(:,2)>point(2)-max_dist&LocMatrixOuttemp(:,1)>point(1)-max_dist);
        testLength(idx)=length(probPoints);
        [indices,dists] = findNearestNeighbors(pointCloud(LocMatrixOuttemp(probPoints,:)),point,k);
%         [~,realidx] = ismember(LocMatrixOuttemp(probPoints(indices),:),LocMatrixOuttemp,'rows'); 
        if ~isempty(indices)
            if length(indices)==k
                idxs(idx*k-k+1:idx*k)=probPoints(indices);
            else
                idxs(idx*k-k+1:idx*k-k+length(indices))=probPoints(indices);
            end 
        end
    end
    idxs(find(idxs==0))=[];
    LocMatrixOuttemp(idxs,:)=[];
    [rows,cols]=size(LocMatrixOuttemp);
    LocMatrixOutClean(i)=mat2cell(LocMatrixOuttemp,[rows],[cols]);
end
toc
%% plot the first image clean outer and innen layer
title_plot='clean bone';
plotTwoImg(cell2mat(LocMatrixOutClean(1)),cell2mat(LocMatrixInClean(1)),title_plot)
legend('outer layer','innen layer')

%% cut off the top part and the bottom part of the bone
% [cutBoneOut] = cutBoneTopBottom(LocMatrixOutClean,z_max_perc,z_min_perc);
% [cutBoneIn] = cutBoneTopBottom(LocMatrixInClean,z_max_perc,z_min_perc);

%% plot the first bone after cutting
% title_plot='cut clean bone';
% plotTwoImg(cell2mat(cutBoneOut(1)),cell2mat(cutBoneIn(1)),title_plot)
% legend('outer layer','innen layer')

%% downsampling and save the downsampled pt images in .ply file
[dSampBoneOut] = dSampPtImg(LocMatrixOutClean,dsamDist_tibia_out);
[dSampBoneIn] = dSampPtImg(LocMatrixInClean,dsamDist_tibia_in);
% save downsampled innen and outer ply file
saveLoc2ply(dSampBoneOut,imgNameSet,dir_path_out,0);
saveLoc2ply(dSampBoneIn,imgNameSet,dir_path_in,1);

%% plot the first bone downsampled
title_plot1='downsampled clean bone outer layer';
title_plot2='downsampled clean bone innen layer';
plotOneImg(cell2mat(dSampBoneOut(1)),title_plot1)
plotOneImg(cell2mat(dSampBoneIn(1)),title_plot2)
