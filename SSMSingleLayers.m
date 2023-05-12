%% prealignment using TMM, cut off around half of the bone, bring the rest half to TMM, then do Principal Component Analysis and show SSM results
% Written by Wuwan Tong 12/05/2023, rwth Aachen

%% set parameters
% Confirm the input images
dirName_r='dSampOut_tibia/r';%'dSampIn_femur/r';%'dSampIn_femur/r';
dirName_l='dSampOut_tibia/l';%'dSampIn_femur/l';
% deside if scale the bone or not in TMM, 1: to scale, 0:not to scale
toScale=0;
% set the number of TMM components and the number of iterate for EM algorithm,and prealignment
M_pre=200;
i_ter_pre=150;
M=3000;
i_ter=300;
% set weight on eigenvector: 3 means 99.865% of the shape, 2 means 97.725%
% of the shape
weightEigVec=3;
% set the weight of standard deviation, e.g. weight=[1 2 3], then plot
% Means-3SD Means-2SD Means-SD Means Means+SD Means+2SD Means+3SD
weightEigVecGrad=1:1:3;
% plot type of SSM visualization, 'points' or 'mesh'
plotType='mesh';
% set the cut off percentage, range:[0,1], if you want to leave 1/3 of the bone, set
% cutOffPerc=2/3
cutOffPerc=0.5;

%% script map
% set parameters
% load images of both left and right side bone, mirror the left bone and rotate 180 grad along z axis, TMMgroupwiseReg does not work well without roration.
% show first 2 bones in dir
% run TMMgroupwiseReg as prealignment
% use transform (from TMM) to do the alignment
% show virtal pts and original pts after alignment
% cut off the top 1/2 (cutOffPerc) of femur
% run TMMgroupwiseReg
% show statistical results  
% SMM visualization on first 3 PC (weight=weightEigVec)
% SMM visualization on first 3 PC, gradient weights (weights=weightEigVecGrad)
%% Brief intorductions of functions
% function [BoneData] = loadImg2Cell(dirpath)
% load all images in dirpath
% input: 
%     dirpath: string
% output:
%     BoneData: 1-by-mun_of_img cell array, each cell stores the location
%     matrix of one point cloud image in 'dirpath'

% function [cutBone] = cut_top(LocMatrixSet,z_max)
% cut the transformed femur point cloud after TMM, cut the z>z_max part
% Input: 
%    LocMatrixSet: 1-by-K cell array, location matrix of point clouds
%    z_max, cut z>z_max part of the femur
% Output:
%    cutBone, 1-by-k cell array of n-by-3 point cloud location, with n<N, the
%    number of rest points after cutting

% function [] = statAnalySSM(SSM)
% show statistical results of SSM

% function [] = plotSSM(SSM,weightEigVec,plotType)
% SMM visualization on each PC
% input: 
%    plotType: 'points' or 'mesh'
%    weightEigVec: weight of standard deviation

% function [] = plotGradSSM(SSM,weightEigVecGrad,plotType)
% gradient weights visualization
% input:
%   plotType: 'points' or 'mesh'
%   weightEigVecGrad: the weight of standard deviation, e.g. weight=[1 2 3], then plot
%                     Means-3SD Means-2SD Means-SD Means Means+SD Means+2SD Means+3SD

%% load images of both left and right side bone, mirror the left bone and rotate 180 grad along z axis, TMMgroupwiseReg does not work well without roration.
% load all images in 'dirName_r' and 'dirName_l' to SyntheticBoneData
BoneData_r = loadImg2Cell(dirName_r);
BoneData_l = loadImg2Cell(dirName_l);
SyntheticBoneData_pre=cell(1,length(BoneData_r)+length(BoneData_l));
SyntheticBoneData_pre(1:length(BoneData_r))=BoneData_r(:);
% SyntheticBoneData_pre(length(BoneData_r)+1:end)=BoneData_l(:);
% mirror the left bone to right shape, and rotate 180 grad alone Z axis
Rmat=[-1 0 0;
      0 -1 0;
      0 0 1];
for i=1:length(BoneData_l)
    mirror_img=cell2mat(BoneData_l(i));
    mirror_img(:,2)=-mirror_img(:,2);
    mirror_img=mirror_img*Rmat;       
    % store the mirrored location matrix of left bone in SyntheticBoneData
    [rows,cols]=size(mirror_img);
    SyntheticBoneData_pre(length(BoneData_r)+i)=mat2cell(mirror_img,[rows],[cols]);
end

%% show first 2 bones in dir
pc1=pointCloud(cell2mat(SyntheticBoneData_pre(1)));
pc2=pointCloud(cell2mat(SyntheticBoneData_pre(2)));
figure
pcshowpair(pc1,pc2)
daspect([1 1 1])
title("SSM-original bones first 2 bone in dir")

%% run TMMgroupwiseReg as prealignment

% downsample to around 1k points to reduce the running time
% dSampSyntheticBoneData_pre = dSampPtImg(SyntheticBoneData_pre,9);

[MU,Transform_pre,TrainingSet_pre,UP,PP,Mcoeffs,nu,convg,SSM_pre]=TMMgroupwiseReg(SyntheticBoneData_pre,M_pre,i_ter_pre,1);

%% test trying to rotate tibia by TMM wothout scaling (tested, failed :(
% SyntheticBoneData_pre1=cell(1,length(SyntheticBoneData_pre));
% for i=1:length(SyntheticBoneData_pre)
%     tempbone=TrainingSet_pre(i).TransfPts;
%     [rows,cols]=size(tempbone);
%     SyntheticBoneData_pre1(i)=mat2cell(tempbone,[rows],[cols]);
% end
% SyntheticBoneData_pre2=cut_top(SyntheticBoneData_pre1,cutOffPerc,toScale);
% [MU,Transform_pre1,TrainingSet_pre1,UP,PP,Mcoeffs,nu,convg,SSM_pre1]=TMMgroupwiseReg_noScale(SyntheticBoneData_pre2,M_pre,i_ter_pre,1);

%% use transform to do the alignment
SyntheticBoneData_temp=cell(1,length(SyntheticBoneData_pre));
centroid = [0,0,0];
for k=1:length(SyntheticBoneData_pre)
    pts = SyntheticBoneData_pre{k};
    cpts = mean(pts,1);
    diffp = cpts-centroid;
    npts = bsxfun(@minus,pts,diffp);
    SyntheticBoneData_temp{k} = npts;
end

figure
for i=1:length(SyntheticBoneData_pre)
    cX = SyntheticBoneData_temp{i};
    rk = Transform_pre(i).R;
    tk = Transform_pre(i).t;
    if toScale
        sk = Transform_pre(i).s;
    elseif ~toScale
        sk =1;
    end
    Xo = bsxfun(@minus,cX,tk);
    tempSR = bsxfun(@times,rk,sk);
    trPts = (tempSR\Xo')'; % Transformed point set
    [rows,cols]=size(trPts);
    SyntheticBoneData_temp(i)=mat2cell(trPts,[rows],[cols]);
    scatter3(trPts(:,1),trPts(:,2),trPts(:,3),'.')
    hold on
end
title('SSM-bones after pre-alignment')
daspect([1 1 1]);

%% show virtal pts and original pts after alignment
pc1=pointCloud(TrainingSet_pre(1).VirtualPts);
pc2=pointCloud(TrainingSet_pre(2).VirtualPts);
figure
pcshowpair(pc1,pc2)
daspect([1 1 1])
title("SSM-first 2 bones in dir after alignment, M=200")
pc1=pointCloud(cell2mat(SyntheticBoneData_temp(1)));
pc2=pointCloud(cell2mat(SyntheticBoneData_temp(2)));
figure
pcshowpair(pc1,pc2)
daspect([1 1 1])
title("SSM-first 2 bones in dir after alignment, original point cloud")
        
%% cut off the top 1/2 of femur

SyntheticBoneData=cut_top(SyntheticBoneData_temp,cutOffPerc,toScale);

%% run TMMgroupwiseReg
if ~toScale
    [MU,Transform,TrainingSet,UP,PP,Mcoeffs,nu,convg,SSM]=TMMgroupwiseReg_noScale(SyntheticBoneData,M,i_ter,1);
elseif toScale
    [MU,Transform,TrainingSet,UP,PP,Mcoeffs,nu,convg,SSM]=TMMgroupwiseReg(SyntheticBoneData,M,i_ter,1);
end
%% show statistical results
statAnalySSM(SSM)

%% SMM visualization on first 3 PC (weight=weightEigVec)
plotSSM(SSM,weightEigVec,plotType);

%% SMM visualization on first 3 PC, gradient weights (weights=weightEigVecGrad)
plotGradSSM(SSM,weightEigVecGrad,plotType);
