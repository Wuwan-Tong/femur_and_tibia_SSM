function [] = plotSSM_compareScale(SSM_scale,SSM_noscale,weightEigVec,plotType)
% SMM visualization on each PC

% input: 
%    plotType: 'points' or 'mesh'
%    weightEigVec: weight of standard deviation


% number of PCs anf Points
numPC=length(SSM_scale.eVals);
numPts=length(SSM_scale.MU)/3;


% scale
%calculate plus and minus on each PC
PCPlus_scale=zeros(3*numPts,numPC);
PCMinus_scale=zeros(3*numPts,numPC);
for i=1:numPC
    PCPlus_scale(:,i)=SSM_scale.MU+(weightEigVec*sqrt(SSM_scale.eVals(i))*SSM_scale.eVecs(:,i))';
    PCMinus_scale(:,i)=SSM_scale.MU-(weightEigVec*sqrt(SSM_scale.eVals(i))*SSM_scale.eVecs(:,i))';
end

% reshape 1-by-(3*numPts) array to numPts-by-3 matrix, for plot
PtsPlus_scale=zeros(numPts,3,numPC);
PtsMinus_scale=zeros(numPts,3,numPC);
PtsMiddle_scale=reshape(SSM_scale.MU,3,[])';
for i=1:numPC
    PtsPlus_scale(:,:,i)=reshape(PCPlus_scale(:,i),3,[])';
    PtsMinus_scale(:,:,i)=reshape(PCMinus_scale(:,i),3,[])';
end

% no scale
%calculate plus and minus on each PC
PCPlus_noscale=zeros(3*numPts,numPC);
PCMinus_noscale=zeros(3*numPts,numPC);
for i=1:numPC
    PCPlus_noscale(:,i)=SSM_noscale.MU+(weightEigVec*sqrt(SSM_noscale.eVals(i))*SSM_noscale.eVecs(:,i))';
    PCMinus_noscale(:,i)=SSM_noscale.MU-(weightEigVec*sqrt(SSM_noscale.eVals(i))*SSM_noscale.eVecs(:,i))';
end

% reshape 1-by-(3*numPts) array to numPts-by-3 matrix, for plot
PtsPlus_noscale=zeros(numPts,3,numPC);
PtsMinus_noscale=zeros(numPts,3,numPC);
PtsMiddle_noscale=reshape(SSM_noscale.MU,3,[])';
for i=1:numPC
    PtsPlus_noscale(:,:,i)=reshape(PCPlus_noscale(:,i),3,[])';
    PtsMinus_noscale(:,:,i)=reshape(PCMinus_noscale(:,i),3,[])';
end


% plot middle, plus and minus for each PCs
% [t,tnorm]=MyRobustCrust(PtsMiddle);
[t] = MyCrustOpen(PtsMiddle_scale); %Triangle connectivity, specified as a 3-column matrix where each row contains the point vertices defining a triangle face.
  
%% plot all pcs in one figure
 
for i=1:3 % plot only first 3 PCs, write for i=1:numPC to plot all PCs 
    figure
    view(0,0)
    % set the title of each subplots
    figNamePCs=append('PC',num2str(i));
    colormap bone
    set(gcf,'Color',[1 1 0.88]) % set color for figure handle   
    % scale
    subplot(2,3,1);
    if strcmp(plotType,'mesh')
        trisurf(t,PtsMiddle_scale(:,1),PtsMiddle_scale(:,2),PtsMiddle_scale(:,3),'Edgecolor','none');
        light
        lighting phong;
        lightangle(0,-90)
    elseif strcmp(plotType,'points')
        scatter3(PtsMiddle_scale(:,1),PtsMiddle_scale(:,2),PtsMiddle_scale(:,3),'.');
    end
    title('mean');
    text(min(xlim)-35,(max(ylim)-min(ylim))*0.6,'Scale')
    set(gca, 'visible', 'off') % axis not wisible
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    axis equal
    
    subplot(2,3,2);
    if strcmp(plotType,'mesh')
        trisurf(t,PtsMinus_scale(:,1,i),PtsMinus_scale(:,2,i),PtsMinus_scale(:,3,i),'Edgecolor','none');
        light
        lighting phong;
        lightangle(0,-90)
    elseif strcmp(plotType,'points')
        scatter3(PtsMinus_scale(:,1,i),PtsMinus_scale(:,2,i),PtsMinus_scale(:,3,i),'.');
    end
    title('-3\sigma');
    text(0,20,max(zlim)+45,figNamePCs);
    set(gca, 'visible', 'off') % axis not wisible
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    axis equal

    subplot(2,3,3);
    if strcmp(plotType,'mesh')
        trisurf(t,PtsPlus_scale(:,1,i),PtsPlus_scale(:,2,i),PtsPlus_scale(:,3,i),'Edgecolor','none');
        light
        lighting phong;
        lightangle(0,-90)
    elseif strcmp(plotType,'points')
        scatter3(PtsPlus_scale(:,1,i),PtsPlus_scale(:,2,i),PtsPlus_scale(:,3,i),'.');
    end
    title('+3\sigma');
    set(gca, 'visible', 'off') % axis not wisible
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    axis equal


    % no scale
    subplot(2,3,4);
    if strcmp(plotType,'mesh')
        trisurf(t,PtsMiddle_scale(:,1),PtsMiddle_scale(:,2),PtsMiddle_scale(:,3),'Edgecolor','none');
        light
        lighting phong;
        lightangle(0,-90)
    elseif strcmp(plotType,'points')
        scatter3(PtsMiddle_scale(:,1),PtsMiddle_scale(:,2),PtsMiddle_scale(:,3),'.');
    end
    % title('mean');
    text(min(xlim)-35,(max(ylim)-min(ylim))*0.6,'Noscale')
    set(gca, 'visible', 'off') % axis not wisible
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    axis equal
    
    subplot(2,3,5);
    if strcmp(plotType,'mesh')
        trisurf(t,PtsMinus_scale(:,1,i),PtsMinus_scale(:,2,i),PtsMinus_scale(:,3,i),'Edgecolor','none');
        light
        lighting phong;
        lightangle(0,-90)
    elseif strcmp(plotType,'points')
        scatter3(PtsMinus_scale(:,1,i),PtsMinus_scale(:,2,i),PtsMinus_scale(:,3,i),'.');
    end
    % title('-3\sigma');
    set(gca, 'visible', 'off') % axis not wisible
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    axis equal

    subplot(2,3,6);
    if strcmp(plotType,'mesh')
        trisurf(t,PtsPlus_scale(:,1,i),PtsPlus_scale(:,2,i),PtsPlus_scale(:,3,i),'Edgecolor','none');
        light
        lighting phong;
        lightangle(0,-90)
    elseif strcmp(plotType,'points')
        scatter3(PtsPlus_scale(:,1,i),PtsPlus_scale(:,2,i),PtsPlus_scale(:,3,i),'.');
    end
    % title('+3\sigma');
    set(gca, 'visible', 'off') % axis not wisible
    set(findall(gca, 'type', 'text'), 'visible', 'on')
    set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[1 1 1]);
    axis equal
end