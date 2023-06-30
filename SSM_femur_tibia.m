%% Apply TMMgroupwiseReg on femur or tibia .ply images, visualize and analyse the resulted SSM
% Written by Wuwan Tong April 2023


%% set parameters
% set downsampling distance
dsamDist=5;
% Confirm the input images
dirName_r='femur_r'; % dirName_r='tibia_r'; for tibia
dirName_l='femur_l';% dirName_l='tibia_l'; for tibia
% set the number of TMM components and the number of iterate for EM algorithm,
M=300;
i_ter=300;
% set weight on eigenvector: 3 means 99.865% of the shape, 2 means 97.725%
% of the shape
weightEigVec=3;
%% load images of both left and right side bone, mirror the left bone and rotate 180 grad along z axis, TMMgroupwiseReg does not work well without roration.
% load all images in 'femur_r' and 'femur_l' or 'tibia_r' and 'tibia_l' to ImgList
ImgList_r=dir(dirName_r);
ImgList_l=dir(dirName_l);

% N: number of images of right bone and left bone
N_r=length(ImgList_r)-2;
N_l=length(ImgList_l)-2;
% create an empty cell array to store point clouds for all femur or tibia 
SyntheticBoneData=cell(1,N_r+N_l);

% i: store image location in i-th cell of the array
i=0;
for n=1:N_r+2
    i=i+1;
    ImgName=ImgList_r(n).name;
    if ImgName=='.'
        i=i-1;
        continue;
    end
    ImgPath=fullfile(dirName_r,ImgName);
    % load images to point cloud and store the location in a matrix
    ptCloud= pcread(ImgPath);
    
    % downsampling to around 2000-3000 points
    DownsampledPtCloud=pcdownsample(ptCloud,'gridAverage',dsamDist);% mode gridAverage, fixed stepsize->different point number
%     pcshow(DownsampledPtCloud);    
    LocMatrix=DownsampledPtCloud.Location;
    
    % store the downsampled matrix(location) in one cell of the array
    [rows,cols]=size(LocMatrix);
    SyntheticBoneData(1,i)=mat2cell(LocMatrix,[rows],[cols]);
end
%deal with femur_l or tibia_l, a mirror flip is required, and roration afterwards
for n=1:N_l+2
    i=i+1;
    ImgName=ImgList_l(n).name;
    if ImgName=='.'
        i=i-1;
        continue;
    end
    ImgPath=fullfile(dirName_l,ImgName);
    % load images to point cloud and store the location in a matrix
    ptCloud= pcread(ImgPath);
    % downsampling to around 2000-3000 points
    DownsampledPtCloud=pcdownsample(ptCloud,'gridAverage',dsamDist);
%     pcshow(DownsampledPtCloud);   
    LocMatrix=DownsampledPtCloud.Location;
    
    % mirror the left bone to right shape, and rotate 180 grad alone Z
    % axis
    Rmat=[-1 0 0;
          0 -1 0;
          0 0 1];
    mirror_img=LocMatrix;
    mirror_img(:,2)=-mirror_img(:,2);
    mirror_img=mirror_img*Rmat;
       
    % store the downsampled matrix(location) in one cell of the array
    [rows,cols]=size(mirror_img);
    SyntheticBoneData(1,i)=mat2cell(mirror_img,[rows],[cols]);
    
end

%% run TMMgroupwiseReg
[MU,Transform,TrainingSet,UP,PP,Mcoeffs,nu,convg,SSM]=TMMgroupwiseReg(SyntheticBoneData,M,i_ter,1);

%% cut the transformed bone and run TMM again
cutBone = cut_femur(TrainingSet,N_r+N_l);
[MU,Transform,TrainingSet,UP,PP,Mcoeffs,nu,convg,SSM]=TMMgroupwiseReg(cutBone,M,i_ter,1);
%% show statistical results
figure
scatter3(SSM.bVecs(:,1),SSM.bVecs(:,2),SSM.bVecs(:,3));
xlabel('component a')
ylabel('component b')
zlabel('component c')
title('position of datapoints on the first 3 PCs')

prop=zeros(1,length(SSM.exp));
for i=1:length(SSM.exp)
    prop(i)=sum(SSM.exp(1:i));
end
figure
yyaxis left
plot(SSM.eVals)
ylabel('eigenvalues')
yyaxis right
plot(prop)
ylabel('percentage')
xlabel('componnents')
xticks(1:1:length(SSM.exp));
ylim([0 100])
title('eigenvalue and percentage distribution of PCs')

%% SMM visualization on each PC
% number of PCs anf Points
numPC=length(SSM.eVals);
numPts=length(SSM.MU)/3;

%calculate plus and minus on each PC
PCPlus=zeros(3*numPts,numPC);
PCMinus=zeros(3*numPts,numPC);
for i=1:numPC
    PCPlus(:,i)=SSM.MU+(weightEigVec*sqrt(SSM.eVals(i))*SSM.eVecs(:,i))';
    PCMinus(:,i)=SSM.MU-(weightEigVec*sqrt(SSM.eVals(i))*SSM.eVecs(:,i))';
end

% reshape 1-by-(3*numPts) array to numPts-by-3 matrix, for plot
PtsPlus=zeros(numPts,3,numPC);
PtsMinus=zeros(numPts,3,numPC);
PtsMiddle=reshape(SSM.MU,3,[])';
for i=1:numPC
    PtsPlus(:,:,i)=reshape(PCPlus(:,i),3,[])';
    PtsMinus(:,:,i)=reshape(PCMinus(:,i),3,[])';
end

% plot middle, plus and minus for each PCs
[t] = MyCrustOpen(PtsMiddle); %Triangle connectivity, specified as a 3-column matrix where each row contains the point vertices defining a triangle face.
for i=1:3 % plot only first 3 PCs, write for i=1:numPC to plot all PCs
    % set the title of each subplots
    figNamePlus=append ('PC',num2str(i),' Plus, (mean+',num2str(weightEigVec),'\sigma)');
    figNameMinus=append ('PC',num2str(i),' Minus,(mean-',num2str(weightEigVec),'\sigma)');   
    figNameMiddle=append ('PC',num2str(i),' Mean');
    
    figure  
    colormap bone
    set(gcf,'Color',[1 1 0.88]) % set color for figure handle   
    view(0,0) % view(az,el) set the viewing angle for a three-dimensional plot, az, is the horizontal rotation 
    % about the z-axis as measured in degrees from the negative y-axis. Positive values indicate counterclockwise 
    % rotation of the viewpoint. el is the vertical elevation of the viewpoint in degrees
%     scatter3(PtsMinus(:,1,i),PtsMinus(:,2,i),PtsMinus(:,3,i),'.');
    subplot(1,3,1);
    trisurf(t,PtsMinus(:,1,i),PtsMinus(:,2,i),PtsMinus(:,3,i),'Edgecolor','none');
    light
    lighting phong;
    set(gca, 'visible', 'off') % axis not wisible
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    axis equal
    title(figNameMinus);

    subplot(1,3,2);
    trisurf(t,PtsMiddle(:,1),PtsMiddle(:,2),PtsMiddle(:,3),'Edgecolor','none');
    light
    lighting phong;
    set(gca, 'visible', 'off') % axis not wisible
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
%     scatter3(PtsMiddle(:,1),PtsMiddle(:,2),PtsMiddle(:,3),'.');
    axis equal
    title(figNameMiddle);

    subplot(1,3,3);
    trisurf(t,PtsPlus(:,1,i),PtsPlus(:,2,i),PtsPlus(:,3,i),'Edgecolor','none');
    light
    lighting phong;
    set(gca, 'visible', 'off') % axis not wisible
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
%     scatter3(PtsPlus(:,1,i),PtsPlus(:,2,i),PtsPlus(:,3,i),'.');
    axis equal
    title(figNamePlus);
end

%% gradient weights plot
% set the weight of standard deviation, e.g. weight=[1 2 3], then plot
% Means-3SD Means-2SD Means-SD Means Means+SD Means+2SD Means+3SD
weightEigVecGrad=1:1:3;

PCPlus=zeros(3*numPts,numPC,length(weightEigVecGrad));
PCMinus=zeros(3*numPts,numPC,length(weightEigVecGrad));
PtsPlus=zeros(numPts,3,numPC,length(weightEigVecGrad));
PtsMinus=zeros(numPts,3,numPC,length(weightEigVecGrad));
for i_weights=1:length(weightEigVecGrad)
    weight=weightEigVecGrad(i_weights);   
    % number of PCs anf Points
    numPC=length(SSM.eVals);
    numPts=length(SSM.MU)/3;
    %calculate plus and minus on each PC, on size 1-by-(3*numPts)
    for i=1:numPC
        PCPlus(:,i,i_weights)=SSM.MU+(weight*sqrt(SSM.eVals(i))*SSM.eVecs(:,i))';
        PCMinus(:,i,i_weights)=SSM.MU-(weight*sqrt(SSM.eVals(i))*SSM.eVecs(:,i))';
    end
    % reshape 1-by-(3*numPts) array to numPts-by-3 matrix, for plot
    PtsMiddle=reshape(SSM.MU,3,[])';
    for i=1:numPC
        PtsPlus(:,:,i,i_weights)=reshape(PCPlus(:,i,i_weights),3,[])';
        PtsMinus(:,:,i,i_weights)=reshape(PCMinus(:,i,i_weights),3,[])';
    end
end
% create cell array to store the titles of subplots
figNamePlusGrad=cell(1,length(weightEigVecGrad));
figNameMinusGrad=cell(1,length(weightEigVecGrad));

[t] = MyCrustOpen(PtsMiddle); %Triangle connectivity, specified as a 3-column matrix where each row contains the point vertices defining a triangle face.
for i=1:3 % plot only first 3 PCs, write for i=1:numPC to plot all PCs
    % set the title of each subplots
    for i_weights=1:length(weightEigVecGrad)
        figNamePlusGrad(i_weights)=cellstr(append ('PC',num2str(i),' +',num2str(weightEigVecGrad(i_weights)),'\sigma'));
        figNameMinusGrad(i_weights)=cellstr(append ('PC',num2str(i),' -',num2str(weightEigVecGrad(i_weights)),'\sigma'));
    end
    figNameMiddle=append ('PC',num2str(i),' Mean');
   
    figure
    colormap bone;
    for i_weights=1:length(weightEigVecGrad)
        % plot minus
        subplot(1,length(weightEigVecGrad)*2+1,i_weights);
        trisurf(t,PtsMinus(:,1,i,length(weightEigVecGrad)-i_weights+1),PtsMinus(:,2,i,length(weightEigVecGrad)-i_weights+1),PtsMinus(:,3,i,length(weightEigVecGrad)-i_weights+1),'Edgecolor','none');
        light
        lighting phong;
        set(gca, 'visible', 'off') % axis not wisible
        set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    %     scatter3(PtsMinus(:,1,i),PtsMinus(:,2,i),PtsMinus(:,3,i),'.');
        axis equal
        title(figNameMinusGrad(length(weightEigVecGrad)-i_weights+1));
        % plot plus
        subplot(1,length(weightEigVecGrad)*2+1,length(weightEigVecGrad)+1+i_weights);
        trisurf(t,PtsPlus(:,1,i,i_weights),PtsPlus(:,2,i,i_weights),PtsPlus(:,3,i,i_weights),'Edgecolor','none');
        light
        lighting phong;
        set(gca, 'visible', 'off') % axis not wisible
        set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    %     scatter3(PtsPlus(:,1,i),PtsPlus(:,2,i),PtsPlus(:,3,i),'.');
        axis equal
        title(figNamePlusGrad(i_weights));
    end   
    % plot mean
    subplot(1,length(weightEigVecGrad)*2+1,length(weightEigVecGrad)+1);
    trisurf(t,PtsMiddle(:,1),PtsMiddle(:,2),PtsMiddle(:,3),'Edgecolor','none');
    light
    lighting phong;
    set(gca, 'visible', 'off') % axis not wisible
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
%     scatter3(PtsMiddle(:,1),PtsMiddle(:,2),PtsMiddle(:,3),'.');
    axis equal
    title(figNameMiddle);
end
  
    



