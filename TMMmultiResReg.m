%% Group-wise multi-resolution rigid registration of point sets using mixture of T's
% Written by Nishant Ravikumar 17/11/2015, The University of Sheffield

% Inputs: 
% Tset {Nk x N} is a cell array containing point sets/shapes to be aligned
% M is the number of mixture components to be used (i.e. number of points
% defining the mean shape)
% max_it is the number of EM iterations
% res is the number of resolution levels to be used
% init is a flag to initialise the registration by first recovering global
% translations and centering the shapes at (0,0,0)

% Outputs: MU --> Mean shape, Transform (struct) --> Similarity transform
% parameters for each shape (R-rotation, t-translation, s-scaling),
% TrainingSet (struct) --> Containing aligned
% point sets stored as TrainingSet.TransfPts and aligned soft
% correspondences stored as TrainingSet.VirtualPts, UP, PP are corrected
% and original posterior probabilities respectively, Mcoeffs --> mixture
% coefficients, nu --> set of degrees of freedom for each t-mixture
% component, convg--> change in mean shape to assess convergence, SSM -->
% Statistical shape model trained by applying PCA over estimated
% correspondences

% Please cite: Ravikumar, N., Gooya, A., Cimen, S., Frangi, A.F. and Taylor, Z.A., 2016, October. 
% A multi-resolution t-mixture model approach to robust group-wise alignment of shapes. 
% In International Conference on Medical Image Computing and Computer-Assisted Intervention (pp. 142-149). 
% Springer International Publishing.

function [MU,Transform,TrainingSet,UP,Mcoeffs,nu,SSM,convg]=TMMmultiResReg(Tset,M,max_it,res,init)
tic
% Number of training shapes
K = size(Tset,2);
data = Tset{1};
[~,D] = size(data);
SSM=[];
figure
for k=1:K
    trPts = Tset{k};
    if D==3
    scatter3(trPts(:,1),trPts(:,2),trPts(:,3),'.')
    hold on
    else
        plot(trPts(:,1),trPts(:,2),'.','MarkerSize',15,'LineWidth',2)
        hold on
    end
end
daspect([1 1 1]);
grid off;
if nargin < 4
    res=3;
end
if nargin < 5
    init=1;
end

if init==1
centroid = [0,0,0];
for k=1:K
    pts = Tset{k};
    cpts = mean(pts,1);
    diffp = cpts-centroid;
    Transform(k).Init = diffp;
    npts = bsxfun(@minus,pts,diffp);
    Tset{k} = npts;
end
end

dof=3;
nu= ones(M,1)*dof;
data = Tset{1};
for k=2:K
    data = cat(1,data,Tset{k});
end
[TN,D] = size(data);

R=repmat(eye(D),[1 1 K]);
t_k=repmat(zeros(1,D),[1 1 K]);
s=repmat(1.0,[1 1 K]);
tempA = eye(D);
Aff = repmat(tempA,[1 1 K]);

[MU,PP,UP,tempU,Mcoeffs,Var2] = InitializeSMM(Tset,data,M,nu,Aff,t_k);
TrainingSet=struct;
if D==3
nD = (TN)*D;
else
    nD = (TN-1)*D;
end
figure
if D==3
H = scatter3(MU(:,1),MU(:,2),MU(:,3),'.');
else
H = plot(MU(:,1),MU(:,2),'.','MarkerSize',15,'LineWidth',2);
end
daspect([1 1 1]);
grid off;
convg=struct;
r=0;
rs=1;
%% Start multi-resolution alignment
while r<res
t=1;
%% EM algorithm to estimate GMM parameters and group-wise reg. parameters (R, s, t)
while t < (max_it+1)
    disp(['EM Iteration = ' num2str(t) ' Resolution level = ' num2str(r+1) ' Number of mean model points = ' num2str(M)])
    clear('PX');
    clear('sPP');
    clear('P1');
    clear('PU');
    % M-step
    for n=1:K
    cPP = PP{n};
    cUP = UP{n};
    A_k = Aff(:,:,n);
    bk = t_k(:,:,n);
    A = MU;
    nX = Tset{n};    
    [cPX,cR,cs,ct_k,csPP] = maximization(nX,cPP,cUP,A);
    PX(:,:,n) = cPX;
    R(:,:,n) = cR;
    s(:,:,n) = cs;
    t_k(:,:,n) = ct_k;
    sPP(:,:,n) = csPP;
    Aff(:,:,n) = cs.*cR;
    Transform(n).R=R(:,:,n); tempT=t_k(:,:,n);tempS=s(:,:,n);    
    Transform(n).s = tempS;
    Transform(n).t=tempT;
    P = UP{n};
    P1(:,:,n) = (cs^2).*sum(P,1)';
    end
    % Compute MU, variance & mixture coefficients
    oldMu=MU; 
    MU=bsxfun(@times,sum(PX,3),1./sum(P1,3));
    dM(t) = norm((oldMu-MU),'fro')./norm(oldMu,'fro');
    if t>30
        id = find((dM(t-4:t)) < 1e-3);
        length(id)
        if length(id)==5
        break;
        end
    end
    disp(['dM = ' num2str(dM(t)) ])
    convg(rs).dM(t) = dM(t);
    ssPP = sum(sPP,3);
    NP=sum(ssPP);
    Mcoeffs = bsxfun(@rdivide,ssPP,NP);
   
    for k=1:K
        nX = Tset{k};
       [N,~] = size(nX);
        A=MU;
        wPP = UP{k};
        tr_k = t_k(:,:,k);
        A_k = Aff(:,:,k);
        b_k = repmat(tr_k,[M 1]);
        tMU = (A_k*A')' + b_k;
        euD = pdist2(nX,tMU,'euclidean');
        Xo = euD.^2;   
        varSum(k) = sum(sum((wPP.*Xo)));
    end
    Var2 = sum(varSum)/(nD)
    
%% Update degrees of freedom (nu)

% Estimate (M) nu parameters
[M,D] = size(MU);
for n=1:K
    U = tempU{n};
    lU = log(U) - U;
    cPP = PP{n};
 PU(:,n) = sum(cPP.*lU,1)';
end

sPU = bsxfun(@rdivide,sum(PU,2),ssPP');
    for i=1:M
    ip=nu(i);
    nu(i) = NewtonRaphson(ip,sPU(i),D);
    end   
 maxb = 100.0;
 minb = 2.1;
 nu=max(min(nu,maxb),minb);
 clear('PP');
 clear('UP');
 clear('tempU');
 
    % E-step
    for m=1:K
    nX = Tset{m};
    tr_k = t_k(:,:,m);
    b_k = repmat(tr_k,[M 1]);
    A_k = Aff(:,:,m);
    A = (A_k*MU')' + b_k;
    [PP{m},UP{m},tempU{m}] = expectation(nX,A,Var2,nu,Mcoeffs);
    end
    
    set(gcf,'color','white')
    axis tight;
    drawnow;
    view(0,90)
    if D==3
    set(H,'XData',MU(:,1),'YData',MU(:,2),'ZData',MU(:,3));
    else
        set(H,'XData',MU(:,1),'YData',MU(:,2));
    end
 
    t=t+1;
end
% Compute virtual correspondences for each training sample based on computed R, t, s
    figure
    daspect([1 1 1]);
    for i=1:K
        wPP = UP{i};
        cX = Tset{i};
        rk = Transform(i).R;
        tk = Transform(i).t;
        sk = Transform(i).s;
        sP = sum(wPP,1)';
        Xo = bsxfun(@minus,cX,tk);
        tempSR = bsxfun(@times,rk,sk);
        trPts = (tempSR\Xo')'; % Transformed point set
        virtPts = bsxfun(@rdivide, (wPP'*trPts), sP); % Compute virtually correspondent points
        TrainingSet(i).VirtualPts = virtPts;
        TrainingSet(i).TransfPts = trPts;
        if D==3
        scatter3(trPts(:,1),trPts(:,2),trPts(:,3),'.')
        hold on
        else
        plot(trPts(:,1),trPts(:,2),'.','MarkerSize',15,'LineWidth',2)
        hold on    
        end
    end
    figure
    for j=1:K
        corPts = TrainingSet(j).VirtualPts;
        if D==3
        scatter3(corPts(:,1),corPts(:,2),corPts(:,3),'.')
        hold on
        else
        plot(corPts(:,1),corPts(:,2),'.','MarkerSize',15,'LineWidth',2)
        hold on    
        end
    end
    daspect([1 1 1]);

%% Adaptive sampling from TMM using multinomial distribution
if r<(res-1)

numS = mnrnd(M,Mcoeffs);
[M,~] = size(MU);
sig = eye(D)*Var2;
newMU = MU;
for j=1:M
    cNU = nu(j);
    if numS(j)~=0
    nP = mvnrnd(zeros(numS(j),D),sig);
    newPts = nP*sqrt(cNU/chi2rnd(cNU));
    nmj = bsxfun(@plus,MU(j,:),newPts); 
    newMU = cat(1,newMU,nmj);
    end
end

clear('MU');
MU=newMU;
[M,~] = size(MU);
disp(['Number of model points = ' num2str(M)])

% Re-initialize parameters
dof=3;
oldNu = nu;
nM = M - length(oldNu);
newNu= ones(nM,1)*dof;
clear('nu');
nu = [oldNu;newNu];
Mcoeffs = ones(1,M) * 1/M;

clear('PP');
clear('UP');
clear('tempU');

for i=1:K
    nX = Tset{i};
    tr_k = t_k(:,:,i);
    b_k = repmat(tr_k,[M 1]);
    A_k = Aff(:,:,i);
    A = (A_k*MU')' + b_k;
    [PP{i},UP{i},tempU{i}] = expectation(nX,A,Var2,nu,Mcoeffs);
end
end
r=r+1;
rs=rs+1;
end
%% Train SSM
for k=1:K
    tmpP = TrainingSet(k).VirtualPts(:,1:3);
    pVec = reshape(tmpP',1,[])';
    pMat(:,k) = pVec;
end
muV = reshape(MU',1,[])';
pMat = cat(2,pMat,muV);
[eVecs, bVecs, eVals, ~, EXPLAINED, pcaMu] = pca(pMat');
SSM.MU = pcaMu;
SSM.eVecs = eVecs;
SSM.bVecs = bVecs;
SSM.eVals = eVals;
SSM.exp = EXPLAINED;
filename = ['mrTMM_Model.mat'];
save(filename,'MU','Mcoeffs','nu','TrainingSet','Transform','SSM','convg')
toc

function [MU,pP,UP,U,Mcoeffs,Var2] = InitializeSMM(Tset,data,M,nu,Aff,t_k)
[TN,~] = size(data);
% [~,MU,sumd,~] = kmeans(data,M,'MaxIter',500);
[~,MU,sumd,~] = kmeans(data,M,'MaxIter',500,'Start','uniform','EmptyAction','singleton');
Var2 = sum(sumd)/(TN)
K = size(Tset,2);
Mcoeffs = ones(1,M) * 1/M;
for n=1:K
    nX = Tset{n};
    tr_k = t_k(:,:,n);
    b_k = repmat(tr_k,[M 1]);
    A_k = Aff(:,:,n);
    A = (A_k*MU')' + b_k;
    [pP{n},UP{n},U{n}] = expectation(nX,A,Var2,nu,Mcoeffs);
end
    
function [PP,UP,U] = expectation(X,MU,Var2,nu,priors)
[N,D] = size(X);
[M,~] = size(MU);

mNu = max(nu);
Var2 = (mNu/(mNu-2))*Var2;
sig2 = eye(D)*Var2;
ed2 = pdist2(X,MU,'mahalanobis',sig2); 
n2 = ed2.^2;
Mh=n2;
nuMh = bsxfun(@plus,Mh,nu');
U = bsxfun(@times, (nu+D)', 1./nuMh); % Covariance weights
nuMh = bsxfun(@times,Mh,1./nu');
repNU = repmat(nu',[N 1]);
qP = (repNU+D)./2;
qdT = log1p(nuMh);
lQ = qP.*qdT;
A = gamma((nu+D).*0.5);
T = gamma(nu.*0.5);
normal = (pi*Var2*nu).^(D/2);
lA = log(A);
lT = log(T);
lN = log(normal);
QT = bsxfun(@plus,lQ,lT');
ATQ = bsxfun(@minus,lA',QT);
logPdf = bsxfun(@minus,ATQ,lN');
lPP = bsxfun(@plus,log(priors),logPdf);
s = logsumexp(lPP,2);
logPP = bsxfun(@minus,lPP,s);
PP=exp(logPP);
logUw = log(U);     % compute log of cov. weights
logUP = bsxfun(@plus, logPP, logUw);
UP = exp(logUP);

%% M-step: Compute rigid registration parameters, R, t, s    
function [PX,R,s,b_k,sPP] = maximization(X,PP,UP,MU)

[~,D] = size(X);
[N,M] = size(PP);
% Update Transformation parameters
% Rotation, scaling and translation
PP_m = sum(UP,2); % Sum over model points
NP=sum(PP_m);
rePP = reshape(UP',1,(N*M));

k1 = kron(X,ones(M,1));
k2 = kron(ones(N,1),MU);
% Compute barycenters
x_k = (1/NP)*rePP*k1;
m_k = (1/NP)*rePP*k2;

Xhat = bsxfun(@minus,X, x_k);
Mhat = bsxfun(@minus,MU, m_k);
k1 = kron(Xhat,ones(M,1));
k2 = kron(ones(N,1),Mhat);
tempK = bsxfun(@times,rePP,k1');
A = tempK*k2;
tempG = bsxfun(@times,rePP,k2');
G = tempG*k2;


% Solve for Rotation, scaling, translation, first compute A-matrix
[W,~,V]=svd(A); 

C=eye(D);
tempDET = det(W)*det(V');
C(end,end)=tempDET; % we need strictly rotation (no reflections)

% Compute rotation
R = double(W*C*V');

numr = trace(A'*R);
denm = trace(G);
% Compute scaling
s = (numr/denm);
Ak = s.*R;
% Compute translation
b_k=x_k-(Ak*m_k')';
Xo = bsxfun(@minus, X, b_k);
% Compute PX, required for updating mu (component means)
RT = s.*R';
RtX = (RT*Xo')'; % Transformed point set

PX = (RtX'*UP)';
sPP = sum(PP,1);