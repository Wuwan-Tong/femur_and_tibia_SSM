%% set params
% input the dir name where the input images are stored
dirNameOut_r='dSampOut_tibia/r';
dirNameOut_l='dSampOut_tibia/l';
dirNameIn_r='dSampIn_tibia/r';
dirNameIn_l='dSampIn_tibia/l';

% set the dir name where to store the clean downsampled innen and outer
% layer, images are stored in 'dir_path_out/r' and 'dir_path_out/l'
dir_path_out_r='dSampOut_tibia_rot/r';
dir_path_in_r='dSampIn_tibia_rot/r';
dir_tform_out_r='tform_tibia_out/r';
dir_tform_in_r='tform_tibia_in/r';
dir_path_out_l='dSampOut_tibia_rot/l';
dir_path_in_l='dSampIn_tibia_rot/l';
dir_tform_out_l='tform_tibia_out/l';
dir_tform_in_l='tform_tibia_in/l';
%% find all images right 
ImgListOut=dir(dirNameOut_r);
ImgListIn=dir(dirNameIn_r);
numImage=length(dirNameOut_r)-2;
i=0;% temp variable to locate the first image in folder
find_fixed_flag=0;
for n=1:numImage+2
    i=i+1;
    ImgName=ImgListOut(n).name;
    if ImgName=='.'
        i=i-1;
        continue;
    end
    ImgPathOut_r=fullfile(dirNameOut_r,ImgName);
    % load images to point cloud and store the location in a matrix
    ptCloudOut= pcread(ImgPathOut_r);
    if i==1&&(find_fixed_flag==0)
        ptCloudOut_fixed=ptCloudOut;
        pcwrite(movingRegOut,fullfile(dir_path_out_r,ImgName));
        find_fixed_flag=1;
        continue;
    end
    [tformOut,movingRegOut] = pcregistericp(ptCloudOut,ptCloudOut_fixed);
    pcwrite(movingRegOut,fullfile(dir_path_out_r,ImgName));
    figure
    pcshowpair(movingRegOut,ptCloudOut_fixed);
    save(fullfile(dir_tform_out_r,strrep(ImgName,'.ply','.mat')),"tformOut");

end
i=0;
find_fixed_flag=0;
for n=1:numImage+2
    i=i+1;
    ImgName=ImgListIn(n).name;
    if ImgName=='.'
        i=i-1;
        continue;
    end
    ImgPathIn_r=fullfile(dirNameIn_r,ImgName);
    % load images to point cloud and store the location in a matrix
    ptCloudIn= pcread(ImgPathIn_r);
    if i==1&&(find_fixed_flag==0)
        ptCloudIn_fixed=ptCloudIn;
        find_fixed_flag=1;
        continue;
    end
    [tformIn,movingRegIn] = pcregistericp(ptCloudIn,ptCloudIn_fixed);
    pcwrite(movingRegIn,fullfile(dir_path_in_r,ImgName));
    save(fullfile(dir_tform_in_r,strrep(ImgName,'.ply','.mat')),"tformIn");
end
%% find all images left
ImgListOut=dir(dirNameOut_l);
ImgListIn=dir(dirNameIn_l);
numImage=length(dirNameOut_l)-2;
i=0;% temp variable to locate the first image in folder
find_fixed_flag=0;
for n=1:numImage+2
    i=i+1;
    ImgName=ImgListOut(n).name;
    if ImgName=='.'
        i=i-1;
        continue;
    end
    ImgPathOut_l=fullfile(dirNameOut_l,ImgName);
    % load images to point cloud and store the location in a matrix
    ptCloudOut= pcread(ImgPathOut_l);
    if i==1&&(find_fixed_flag==0)
        ptCloudOut_fixed=ptCloudOut;
        find_fixed_flag=1;
        continue;
    end
    [tformOut,movingRegOut] = pcregistericp(ptCloudOut,ptCloudOut_fixed);
    pcwrite(movingRegOut,fullfile(dir_path_out_l,ImgName));
    save(fullfile(dir_tform_out_l,strrep(ImgName,'.ply','.mat')),"tformOut");

end
i=0;
find_fixed_flag=0;
for n=1:numImage+2
    i=i+1;
    ImgName=ImgListIn(n).name;
    if ImgName=='.'
        i=i-1;
        continue;
    end
    ImgPathIn_l=fullfile(dirNameIn_l,ImgName);
    % load images to point cloud and store the location in a matrix
    ptCloudIn= pcread(ImgPathIn_l);
    if i==1&&(find_fixed_flag==0)
        ptCloudIn_fixed=ptCloudIn;
        find_fixed_flag=1;
        continue;
    end
    [tformIn,movingRegIn] = pcregistericp(ptCloudIn,ptCloudIn_fixed);
    pcwrite(movingRegIn,fullfile(dir_path_in_l,ImgName));
    save(fullfile(dir_tform_in_l,strrep(ImgName,'.ply','.mat')),"tformIn");
end