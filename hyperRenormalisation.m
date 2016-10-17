function data = hyperRenormalisation (X, F0, maxIters, convThresh)
%% hyperRenormalisation, Ellipse Fitting by Hyper-Renormalisation
%
% data = hyperRenormalisation (X)
% data = hyperRenormalisation (X, F0, maxIters, convThresh)
%
% Details This function contains an implementation of the
%     ellipse fitting algorithm as originally described in
%     Hyper-renormalization: Non-minimization approach for
%     geometric estimation by Kanatani, Al-Sharadqah, Chernov
%     and Sugaya (2014). Available here:
%     http://dx.doi.org/10.2197/ipsjtcva.6.143.
%
% Dependencies  This function requires five small functions for
%     calculating the geomtric translation and rotation of the
%     ellipse from the algebriac parameters: "imconicdiscr",
%     "imconicrotate", "imconicrotation", "imconictranslate" and
%     "imconictranslation" from "Fitting quadratic curves and
%     surfaces" by Levente Hunyadi (2014) available from:
%     http://uk.mathworks.com/matlabcentral/fileexchange/45356-fitting-quadratic-curves-and-surfaces
%
% Inputs  X - an n-by-2 array [x,y] of points in 2D
%           F0 - a fixed constant as used in Eq. 3 of the original
%             paper; this can be appropximated as the maximum
%             expected size
%           maxIters - maximum number of iterations to allow
%             (see also convThresh)
%           convThresh - convergence threshold at which to stop
%             iterating (see also maxIters)
%
% Outputs data - a matrix of the ellipse measurements that contains
%               five elements: the center of the ellipse, its major and
%               minor axes and orientation of the major axis.
%
% Examples:
% data = hyperRenormalisation([x,y],100,500,0.01), runs ellipse
%   fitting by hyper-renormalization on the points described by
%   x and y, with an expected size of 100 pixels and stopping
%   iterations when the convergence measure drops below 0.01
%   or the number of iterations exceeds 500.
%
% Copyright 2016 Carl J. Nelson, Durham University, UK
% The original concepts and ideas remain the property of the
% authors: Kanatani, Al-Sharadqah, Chernov and Sugaya.
%
% License   See included LICENSE file or visit the GitHub
%   repository at https://github.com/ChasNelson1990/Hyper-renormalisation-for-Ellipse-Fitting

%% Inputs
if nargin<4 || isempty(convThresh); convThresh = 0.1; end
if nargin<3 || isempty(maxIters); maxIters = 100; end
if nargin<2 || isempty(F0); F0 = max(sqrt(X(:,1).^2 + X(:,2).^2)); end

%% Normalize Data
mx = mean(X(:,1));
my = mean(X(:,2));
sx = (max(X(:,1))-min(X(:,1)))/2;
sy = (max(X(:,2))-min(X(:,2)))/2;

sx = max(sx,sy);
sy = max(sx,sy);

X(:,1) = (X(:,1)-mx)/sx;
X(:,2) = (X(:,2)-my)/sy;

%% Set-Up
npoints = size(X,1);
F0 = F0*ones(npoints,1);

%% Calculate XI
XI = [X(:,1).^2, 2*X(:,1).*X(:,2), X(:,2).^2, 2*X(:,1).*F0, 2*X(:,2).*F0, F0.^2];

%% Calculate V0
V0 = [   4*X(:,1).^2,       4*X(:,1).*X(:,2),           zeros(npoints,1),   4*X(:,1).*F0,       zeros(npoints,1),   zeros(npoints,1),...
         4*X(:,1).*X(:,2),  4*(X(:,1).^2 + X(:,2).^2),  4*X(:,1).*X(:,2),   4*X(:,2).*F0,       4*X(:,1).*F0,       zeros(npoints,1),...
         zeros(npoints,1),  4*X(:,1).*X(:,2),           4*X(:,2).^2,        zeros(npoints,1),	4*X(:,2).*F0,       zeros(npoints,1),...
         4*X(:,1).*F0,      4*X(:,2).*F0,               zeros(npoints,1),   4*F0.^2,            zeros(npoints,1),   zeros(npoints,1),...
         zeros(npoints,1),  4*X(:,1).*F0,               4*X(:,2).*F0,       zeros(npoints,1),   4*F0.^2,            zeros(npoints,1),...
         zeros(npoints,1),  zeros(npoints,1),           zeros(npoints,1),	zeros(npoints,1),   zeros(npoints,1),   zeros(npoints,1)];
V0 = reshape(V0,npoints,6,6);

%% Calculate W
W = ones(npoints,1);

%% Iterations
iters=0;
finish=false;
uold = [0,0,0,0,0,0]';
while ~finish || iters<maxIters
    iters = iters+1;

    M = zeros(6,6);
    for n=1:npoints
        M = M + (W(n) * (XI(n,:)' * XI(n,:)));
    end

    M = M/npoints;

    N1 = zeros(6,6);
    N2 = zeros(6,6);
    e = [ 1 ; 0 ; 1 ; 0 ; 0 ; 0 ]';

    [evec,ev] = eig(M);
    Minv = zeros(6,6);
    for n=2:6
        Minv = Minv + ((1/ev(n,n)) * (evec(:,n)*evec(:,n)'));
    end

    for n=1:npoints
        work = XI(n,:)' * e;
        N1 = N1 + (W(n) * (squeeze(V0(n,:,:)) + (work + work')));
        work = squeeze(V0(n,:,:)) * Minv * (XI(n,:)' * XI(n,:));
        N2 = N2 + (W(n)^2 * (dot(XI(n,:)',(Minv * XI(n,:)')) * squeeze(V0(n,:,:)) + (work + work')));
    end
    N = N1/npoints - N2/(npoints^2);

    [~,~,nvec] = eig(N,M);
    u = nvec(:,6);
    if u(1)<0; u=-u; end;
    u = u/(sqrt(sum(u.^2)));

    if norm(u-uold)<convThresh
        finish = true;
    else
        for n=1:npoints
            W(n) = 1/dot(u,squeeze(V0(n,:,:))*u);
        end
            uold = u;
    end
end
u = u';

u = [u(1); 2*u(2); u(3); 2*u(4); 2*u(5); u(6)];

%% Unnormalise
par = [ u(1)*sy*sy,   ...
        u(2)*sx*sy,   ...
        u(3)*sx*sx,   ...
        -2*u(1)*sy*sy*mx - u(2)*sx*sy*my + u(4)*sx*sy*sy,   ...
        -u(2)*sx*sy*mx - 2*u(3)*sx*sx*my + u(5)*sx*sx*sy,   ...
        u(1)*sy*sy*mx*mx + u(2)*sx*sy*mx*my + u(3)*sx*sx*my*my - u(4)*sx*sy*sy*mx - u(5)*sx*sx*sy*my + u(6)*sx*sx*sy*sy]';

%% Extract Ellipse Parameters
% Adjust parameters in vector p
par(2) = 0.5 * par(2);
par(4:5) = 0.5 * par(4:5);

% Use implicit equation to compute ellipse semi-axes
q = 2*(par(1)*par(5)^2+par(3)*par(4)^2+par(6)*par(2)^2-2*par(2)*par(4)*par(5)-par(1)*par(3)*par(6))/(par(2)^2-par(1)*par(3));
r = realsqrt((par(1)-par(3))^2+4*par(2)^2);
data(:,3) = 2*F0(1)*real(sqrt(q/(r-(par(1)+par(3)))));   % major axis
data(:,4) = 2*F0(1)*real(sqrt(q/(-r-(par(1)+par(3)))));  % minor axis

% Check
if data(:,3)<data(:,4)
    data(:,3) = data(:,4); %major axis
    data(:,4) = 2*F0(1)*real(sqrt(q/(r-(par(1)+par(3)))));%minor axis
end

% Readjust parameters in vector p
par(2) = 2 * par(2);
par(4:5) = 2 * par(4:5);

% Calculate translation and rotation of ellipse
data(:,1:2) = imconictranslation(par);%position in space
R = imconicrotation(imconictranslate(par, -data(:,1:2)/2));%rotation matrix
data(:,5) = acosd(-R(1,1));%orientation of major axis compared to x axis

end
