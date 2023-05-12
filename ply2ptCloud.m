% convert ply images to point cloud and store the location in each cell of
% a cell array SyntheticFemurData

% set downsampleing distance
dsamDist=5;%1.3;
% load all images in 'femur_r' and 'femur_l' to ImgList
ImgList_r=dir('femur_r');
ImgList_l=dir('femur_l');
% N: number of images of femur right and left
N_r=length(ImgList_r)-2;
N_l=length(ImgList_l)-2;

% create an empty cell array to store point clouds for all femur right
SyntheticFemurData=cell(1,N_r+N_l);
% i: store image location in i-th cell of the array
i=0;
for n=1:N_r+2
    i=i+1;
    ImgName=ImgList_r(n).name;
    if ImgName=='.'
        i=i-1;
        continue;
    end
    ImgPath=fullfile('femur_r',ImgName);
    % load images to point cloud and store the location in a matrix
    ptCloud= pcread(ImgPath);
    
    % downsampling to around 1157 points
    % test mode random->got bad result
%     DownsampledPtCloud=pcdownsample(ptCloud,'random',0.014); 
    % mode gridAverage, fixed stepsize->different point number
    DownsampledPtCloud=pcdownsample(ptCloud,'gridAverage',dsamDist);
%     pcshow(DownsampledPtCloud);
    
    LocMatrix=DownsampledPtCloud.Location;
    % store the downsampled matrix(location) in one cell of the array
    [rows,cols]=size(LocMatrix);
    SyntheticFemurData(1,i)=mat2cell(LocMatrix,[rows],[cols]);
end
%% deal with femur_l, a mirror flip is required
for n=1:N_l+2
    i=i+1;
    ImgName=ImgList_l(n).name;
    if ImgName=='.'
        i=i-1;
        continue;
    end
    ImgPath=fullfile('femur_l',ImgName);
    % load images to point cloud and store the location in a matrix
    ptCloud= pcread(ImgPath);

    % downsampling to around 1157 points
    % test mode random->got bad result
%     DownsampledPtCloud=pcdownsample(ptCloud,'random',0.014); 
    % mode gridAverage, fixed stepsize->different point number
    DownsampledPtCloud=pcdownsample(ptCloud,'gridAverage',dsamDist);
%     pcshow(DownsampledPtCloud);
    
    LocMatrix=DownsampledPtCloud.Location;
    
    % mirror the left femur to right shape, and rotate 180 grad
    D=[-1 0 0;
       0 -1 0;
       0 0 1];
    mirror_img=LocMatrix;
    mirror_img(:,2)=-mirror_img(:,2);
    mirror_img=mirror_img*D;
       
    % store the downsampled matrix(location) in one cell of the array
    [rows,cols]=size(mirror_img);
    SyntheticFemurData(1,i)=mat2cell(mirror_img,[rows],[cols]);
    
%     % test start-----
%     mirror_img=LocMatrix;
%     mirror_img(:,2)=-mirror_img(:,2);
%     centroid = [0,0,0];
% 
%     cmirror_img = mean(mirror_img,1);
%     diffp = cmirror_img-centroid;
%     npts = bsxfun(@minus,mirror_img,diffp);
%     mirror_img = npts;
% 
%     figure
%     scatter3(mirror_img(:,1),mirror_img(:,2),mirror_img(:,3))
%     figure
%     scatter3(LocMatrix(:,1),LocMatrix(:,2),LocMatrix(:,3))
%     %test end-----

end

