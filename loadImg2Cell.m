function [BoneData] = loadImg2Cell(dirpath)
% load all images in dirpath
ImgList=dir(dirpath);
% N: number of images of right bone and left bone
N=length(ImgList)-2;
% create an empty cell array to store point clouds for all femur or tibia 
BoneData=cell(1,N);
% i: store image location in i-th cell of the array
i=0;
for n=1:N+2
    i=i+1;
    ImgName=ImgList(n).name;
    if ImgName=='.'
        i=i-1;
        continue;
    end
    ImgPath=fullfile(dirpath,ImgName);
    % load images to point cloud and store the location in a matrix
    ptCloud= pcread(ImgPath);      
    LocMatrix=ptCloud.Location;
    % store the downsampled matrix(location) in one cell of the array
    [rows,cols]=size(LocMatrix);
    BoneData(1,i)=mat2cell(LocMatrix,[rows],[cols]);
end
end
