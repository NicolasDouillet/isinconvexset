function isin = isinconvexset(V, H, M)
%% isinconvexset : function to check if a vertex is located inside or outside a given
% convex set, boundary not included (opened set). Supports dimensions 2 and 3.
%
% Author & support : nicolas.douillet (at) free.fr, 2018-2023.
%
%
% Syntax
%
% isin = isinconvexset(V, H, M);
%
%
% Description
%
% isin = isinconvexset(V, H, M) computes the boolean isin which is true/logical 1
% in the case the vertex M belongs to the opened convex set (V,H) defined by
% the vertex set V, and its hyperplane set H. isin is false/logical 0 in the
% case vertex M belongs to the complementary set or the boundary convex hull.
%
%
% See also : CONVHULL
%
%
% Input arguments
%
%       [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the convex set, with size(V,1) > size(V,2) to define a relevant convex set.
%       [ |  |  |]
%
%       [ |  |]
% - H = [i1 i2], positive integer matrix double, the hyperplane set (edge list in 2D, face list in 3D). 
%       [ |  |]
%
%                H elements / rows don't need to be coherently oriented. Size(H) = [nb_hyperplanes,2].
%
%       [ |  |  |]
% - M = [Mx My Mz], real row vector or matrix double, the coordinates of the vertex / vertices  to check. Size(M,2) = size(V,2).
%       [ |  |  |]
%
%
% Output argument
%
%          [      |      ]
% - isin = [logical 1 / 0], logical true (1)/false (0) scalar / column vector. The boolean result. Size(isin) = [size(M,1),1].
%          [      |      ]
%
%
% Example #1 : 2D convex hull of random point cloud
% N = 28;
% V = 2*(rand(N,2)-0.5);
%
% H_raw = convhull(V);
% nh = numel(H_raw);
% H = repelem(H_raw',cat(2,1,2*ones(1,nh-2),1));
% H = reshape(H,[2,nh-1])';
%
% [A B] = meshgrid(-1:0.1:1);
% M = cat(2,A(:),B(:));
%
% isin = isinconvexset(V,H,M);
%
% figure
% line(V(H_raw,1),V(H_raw,2),'Color',[1 1 0],'LineWidth',2), hold on;
% ColorSpec = num2cell(cat(2,~isin,isin,zeros(size(isin,1),1)),2);
% cellfun(@(r1,r2) plot(r1(1,1),r1(1,2),'+','Color',r2,'MarkerSize',4,'LineWidth',4),num2cell(M,2),ColorSpec,'un',0);
% set(gcf,'Color',[0 0 0]), set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'FontSize',16);
% xlabel('X'), ylabel('Y');
% axis equal, axis tight;
%
%
% Example #2 : 3D convex hull of random point cloud
% N = 32;
% V = 2*(rand(N,3)-0.5);
%
% [A B C] = meshgrid(-1:0.3:1,-1:0.3:1,-1:0.3:1);
% M = cat(2,A(:),B(:),C(:));
%
% H = convhull(V);
% vtx_idx = unique(H(:));
%
% isin = isinconvexset(V,H,M);
%
% figure;
% trisurf(H,V(:,1),V(:,2),V(:,3)), shading faceted, hold on;
% colormap([0 1 1]);
% ColorSpec = num2cell(cat(2,~isin,isin,zeros(size(isin,1),1)),2);
% cellfun(@(r1,r2) plot3(r1(1,1),r1(1,2),r1(1,3),'+','Color',r2,'MarkerSize',4,'LineWidth',4),num2cell(M,2),ColorSpec,'un',0);
% set(gcf,'Color',[0 0 0]), set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',[1 1 1],'FontSize',16);
% ax = gca;
% ax.Clipping = 'off';
% alpha(0.4);
% xlabel('X'), ylabel('Y'), zlabel('Z');
% axis equal, axis tight;
% camlight left;


%% Input parsing
assert(nargin > 2,'Not enought input arguments.');
assert(nargin < 4,'Too many input arguments.');
assert(isequal(size(V,2),size(H,2),size(M,2)),'All the inputs must have the same number of colums (two or three here).');


%% Body
% epsilon = 1e4*eps;
[~,dim] = size(V);
G = mean(V,1);
nb_h = size(H,1);

% Hyperplane isobarycentres
Gh = cell2mat(cellfun(@(r) mean(V(r,:),1),num2cell(H,2),'un',0));

if dim == 2 % 2D
    
    % Face normals
    H1 = V(H(:,1),:); % X
    H2 = V(H(:,2),:); % Y
    n = [H1(:,2)-H2(:,2), H2(:,1)-H1(:,1)]; % = [yA-yB; xB-xA]        
    
elseif dim == 3 % 3D
    
    % Face normals
    n = cross(V(H(:,2),:)-V(H(:,1),:),V(H(:,3),:)-V(H(:,1),:));
    
    % False 3D case detection    
    
else
    
    error('isinconvexset only works for convex sets of dimension 2 or 3.');
    
end

% Vectors from hyperplane barycentres to candidate points
Gv = cellfun(@(r) repmat(r,[size(Gh,1),1]) - Gh,num2cell(M,2),'un',0);

% Orient coherently
orientation = sign(dot(n,repmat(G,[nb_h,1])-Gh,2));

if ~isequal(orientation,-ones(nb_h,1)) && ~isequal(orientation,ones(nb_h,1))
    
    n = n.*orientation;
    
end

% Signed dot product
sdot_prod = cellfun(@(c) sign(dot(c,n,2)),Gv,'un',0);

% Compute isin logical value
isin = cell2mat(cellfun(@(c) isequal(c,ones(nb_h,1)) || ...
                             isequal(c,-ones(nb_h,1)),sdot_prod,'un',0));


end % isinconvexset