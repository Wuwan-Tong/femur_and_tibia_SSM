function [] = plotGradSSM(SSM,weightEigVecGrad,plotType)
% gradient weights visualization

% input:
%   plotType: 'points' or 'mesh'
%   weightEigVecGrad: the weight of standard deviation, e.g. weight=[1 2 3], then plot
%                     Means-3SD Means-2SD Means-SD Means Means+SD Means+2SD Means+3SD


% number of PCs anf Points
numPC=length(SSM.eVals);
numPts=length(SSM.MU)/3;

PCPlus=zeros(3*numPts,numPC,length(weightEigVecGrad));
PCMinus=zeros(3*numPts,numPC,length(weightEigVecGrad));
PtsPlus=zeros(numPts,3,numPC,length(weightEigVecGrad));
PtsMinus=zeros(numPts,3,numPC,length(weightEigVecGrad));
for i_weights=1:length(weightEigVecGrad)
    weight=weightEigVecGrad(i_weights);   
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
    figNameMiddle=append ('Mean');
   
    figure
    colormap bone;
    for i_weights=1:length(weightEigVecGrad)
        % plot minus
        subplot(1,length(weightEigVecGrad)*2+1,i_weights);
        if strcmp(plotType,'mesh')
            trisurf(t,PtsMinus(:,1,i,length(weightEigVecGrad)-i_weights+1),PtsMinus(:,2,i,length(weightEigVecGrad)-i_weights+1),PtsMinus(:,3,i,length(weightEigVecGrad)-i_weights+1),'Edgecolor','none');
            light
            lighting phong;
        elseif strcmp(plotType,'points')
            scatter3(PtsMinus(:,1,i),PtsMinus(:,2,i),PtsMinus(:,3,i),'.');
        end
        title(figNameMinusGrad(length(weightEigVecGrad)-i_weights+1));
        set(gca, 'visible', 'off') % axis not wisible
        set(findall(gca, 'type', 'text'), 'visible', 'on')
        set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);   
        axis equal
        
        % plot plus
        subplot(1,length(weightEigVecGrad)*2+1,length(weightEigVecGrad)+1+i_weights);
        if strcmp(plotType,'mesh')
            trisurf(t,PtsPlus(:,1,i,i_weights),PtsPlus(:,2,i,i_weights),PtsPlus(:,3,i,i_weights),'Edgecolor','none');
            light
            lighting phong;
        elseif strcmp(plotType,'points')
            scatter3(PtsPlus(:,1,i),PtsPlus(:,2,i),PtsPlus(:,3,i),'.');
        end
        title(figNamePlusGrad(i_weights));
        set(gca, 'visible', 'off') % axis not wisible
        set(findall(gca, 'type', 'text'), 'visible', 'on')
        set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]); 
        axis equal
        
    end   
    % plot mean
    subplot(1,length(weightEigVecGrad)*2+1,length(weightEigVecGrad)+1);
    if strcmp(plotType,'mesh')
        trisurf(t,PtsMiddle(:,1),PtsMiddle(:,2),PtsMiddle(:,3),'Edgecolor','none');
        light
        lighting phong;
     elseif strcmp(plotType,'points')
         scatter3(PtsMiddle(:,1),PtsMiddle(:,2),PtsMiddle(:,3),'.');
    end
    title(figNameMiddle);
    set(gca, 'visible', 'off') % axis not wisible
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    axis equal
    
end


end