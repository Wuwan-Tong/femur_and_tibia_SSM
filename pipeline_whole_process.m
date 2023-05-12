%% pipeline to do SSM from the raw dataset we have
% Written by Wuwan Tong 12/05/2023, rwth Aachen

%% README!!!:
% If you use this pipeline script to run the whole process, please comment
% out the first sections (set parameters) in both file 'cut_innen.m' and
% 'SSMSingleLayers.m', and set parameters in this script instead.

%% required structure of folders:
% to store the input .ply images:
% -aussen_femur
% -innen_femur
% -aussen_tibia
% -innen_tibia

% to store the interim .ply images (clean downsampled single layer bone after
% cut), these images will be used for SSM. 
% you only need to create empty folders, the images will be saved automatically by the code
% -dSampOut_femur
% ---r
% ---l
% -dSampIn_femur
% ---r
% ---l
% -dSampOut_tibia
% ---r
% ---l
% -dSampIn_tibia
% ---r
% ---l

% please make sure these folders are in the same folder as the code,
% otherwise please give the full path in the code

%% map of this script, you can run each section individually

% femur--devide innen and outer layer

% femur--set universal parameters for SSM
% femur--SSM, outer layer, with scale
% femur--SSM, outer layer, without scale
% femur--SSM, innen layer, with scale
% femur--SSM, innen layer, without scale

% tibia--devide innen and outer layer

% tibia--set universal parameters for SSM
% tibia--SSM, outer layer, with scale
% tibia--SSM, outer layer, without scale
% tibia--SSM, innen layer, with scale
% tibia--SSM, innen layer, without scale

%% femur, devide innen and outer layer 
% input the dir name where the input images are stored
dirNameTwolayer='aussen_femur';
dirNameIn='innen_femur';
% set downsampling distance
dsamDist_femur_out=3.4;
dsamDist_femur_in=2.4;
% set the dir name where to store the clean downsampled innen and outer
% layer, images are stored in 'dir_path_out/r' and 'dir_path_out/l'
dir_path_out='dSampOut_femur';
dir_path_in='dSampIn_femur';
% cut off the top and bottom part, for now we do cutting after prealignment instead, so the 
% code for cutting in this part is commented out, you can reuse this part
% later if it's required
% z_max_perc=0.5; % range:[0,1], the region z>z_max will be cut off, 0.6 means top 40% of the bone will be cut off
% z_min_perc=0.1; % range:[0,1], the region z<z_min will be cut off, 0.1 means bottom 10% of the bone will be cut off

% RUN
run("cut_innen.m")

%% femur, do SSM for 4 scenarios: 1.outer layer, with scale; 2.outer layer, without scale; 3.innen layer, with scale; 4.innen layer, without scale
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
% set the cut off percentage, if you want to leave 1/3 of the bone, set
cutOffPerc=0.5;

%% femur, 1.outer layer, with scale
% Confirm the input images
dirName_r='dSampOut_femur/r';%'dSampIn_femur/r';%'dSampIn_femur/r';
dirName_l='dSampOut_femur/l';%'dSampIn_femur/l';
% deside if scale the bone or not in TMM, 1: to scale, 0:not to scale
toScale=1;

% RUN
run("SSMSingleLayers.m")

%% femur, 2.outer layer, without scale
% Confirm the input images
dirName_r='dSampOut_femur/r';
dirName_l='dSampOut_femur/l';
% deside if scale the bone or not in TMM, 1: to scale, 0:not to scale
toScale=0;

% RUN
run("SSMSingleLayers.m")

%% femur, 3.innen layer, with scale
% Confirm the input images
dirName_r='dSampIn_femur/r';
dirName_l='dSampIn_femur/l';
% deside if scale the bone or not in TMM, 1: to scale, 0:not to scale
toScale=1;

% RUN
run("SSMSingleLayers.m")

%% femur, 4.innen layer, without scale
% Confirm the input images
dirName_r='dSampIn_femur/r';
dirName_l='dSampIn_femur/l';
% deside if scale the bone or not in TMM, 1: to scale, 0:not to scale
toScale=0;

% RUN
run("SSMSingleLayers.m")


%% tibia, devide innen and outer layer 
% input the dir name where the input images are stored
dirNameTwolayer='aussen_tibia';
dirNameIn='innen_tibia';
% set downsampling distance
dsamDist_tibia_out=2.5;
dsamDist_tibia_in=1.7;
% set the dir name where to store the clean downsampled innen and outer
% layer, images are stored in 'dir_path_out/r' and 'dir_path_out/l'
dir_path_out='dSampOut_tibia';
dir_path_in='dSampIn_tibia';
% cut up and bottom part, for now we do cutting after prealignment, so the 
% code for cutting in this part is commented out, you can use this part
% later if it's required
z_max_perc=0.5; % the region z>z_max will be cut off, 0.6 means top 40% of the bone will be cut off
z_min_perc=0.1; % the region z<z_min will be cut off, 0.1 means bottom 10% of the bone will be cut off

% RUN
run("cut_innen.m")

%% tibia, do SSM for 4 scenarios: 1.outer layer, with scale; 2.outer layer, without scale; 3.innen layer, with scale; 4.innen layer, without scale
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
% set the cut off percentage, if you want to leave 1/3 of the bone, set
cutOffPerc=0.5;

%% tibia, 1.outer layer, with scale
% Confirm the input images
dirName_r='dSampOut_tibia/r';
dirName_l='dSampOut_tibia/l';
% deside if scale the bone or not in TMM, 1: to scale, 0:not to scale
toScale=1;

% RUN
run("SSMSingleLayers.m")

%% tibia, 2.outer layer, without scale
% Confirm the input images
dirName_r='dSampOut_tibia/r';
dirName_l='dSampOut_tibia/l';
% deside if scale the bone or not in TMM, 1: to scale, 0:not to scale
toScale=0;

% RUN
run("SSMSingleLayers.m")

%% tibia, 3.innen layer, with scale
% Confirm the input images
dirName_r='dSampIn_tibia/r';
dirName_l='dSampIn_tibia/l';
% deside if scale the bone or not in TMM, 1: to scale, 0:not to scale
toScale=1;

% RUN
run("SSMSingleLayers.m")

%% femur, 4.innen layer, without scale
% Confirm the input images
dirName_r='dSampIn_tibia/r';
dirName_l='dSampIn_tibia/l';
% deside if scale the bone or not in TMM, 1: to scale, 0:not to scale
toScale=0;

% RUN
run("SSMSingleLayers.m")

