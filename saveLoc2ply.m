function [] = saveLoc2ply(LocMatrixSet,imgNameSet,dir_path,isInnen)
% find the number of images
numImage=length(LocMatrixSet);
% convert imgNameSet, the cell array to string array
imgNameStr=char(imgNameSet);

for i=1:numImage
    % add innen or aussen to filename
    if isInnen
        newImgName=strrep(imgNameStr(i,:),'.','_innen.');
    else
        newImgName=strrep(imgNameStr(i,:),'.','_aussen.');
    end
    % save pt images
    ptImg=pointCloud(cell2mat(LocMatrixSet(i)));
    if contains(newImgName,'_r') 
        filename=strtrim(fullfile(dir_path,'r',newImgName));
    else
        filename=strtrim(fullfile(dir_path,'l',newImgName));
    end
    pcwrite(ptImg,filename)
end
end