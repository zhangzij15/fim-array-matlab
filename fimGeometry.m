function varargout = fimGeometry(action,varargin)
%FIMGEOMETRY Fixed reference, rigid subarrays and in-plane movable phase centres.
switch action
    case 'prepare'
        cfg = varargin{1};
        if strcmp(cfg.deform.family,'subarray')
            varargout{1} = sa_prepareModel(cfg);
        elseif strcmp(cfg.deform.family,'pca')
            varargout{1} = pc_prepareModel(cfg);
        else
            error('FIM:UnsupportedGeometry','Only subarray and pca are supported.');
        end
    case 'orientation'
        varargout{1}=buildOrientationGeometry(varargin{:});
    case 'build'
        model = varargin{1};
        if strcmp(model.family,'subarray')
            [varargout{1},varargout{2},varargout{3}] = sa_buildGeometry(varargin{:});
        elseif strcmp(model.family,'pca')
            [varargout{1},varargout{2},varargout{3}] = pc_buildGeometry(varargin{:});
        else
            error('FIM:UnsupportedGeometry','Only subarray and pca are supported.');
        end
    otherwise
        error('Unknown fimGeometry action: %s',action);
end
end

function model = sa_prepareModel(cfg)
assert(isfield(cfg,'subarray') && isstruct(cfg.subarray),'Missing cfg.subarray.');
s = cfg.subarray;
sa_checkPositiveInteger(cfg.Nh,'Nh');
sa_checkPositiveInteger(cfg.Nv,'Nv');
sa_checkPositiveInteger(s.m,'subarray.m');
sa_checkPositiveInteger(s.n,'subarray.n');
assert(mod(cfg.Nh,s.m)==0 && mod(cfg.Nv,s.n)==0,'subarray.m/n must divide Nh/Nv respectively; unequal edge blocks are not supported.');
assert(isnumeric(s.min_spacing_lambda) && isreal(s.min_spacing_lambda) && isscalar(s.min_spacing_lambda) && isfinite(s.min_spacing_lambda) && s.min_spacing_lambda>=0,'Invalid subarray.min_spacing_lambda.');
grids = {sa_withZero(s.displacement_grid_lambda),sa_withZero(s.yaw_grid_deg),sa_withZero(s.tilt_grid_deg)};
assert(all(abs(grids{2})<=180) && all(abs(grids{3})<=90),'Subarray yaw/tilt ranges must lie within +/-180 and +/-90 degrees.');
lambda = cfg.c0/cfg.fc_Hz;
x = ((0:cfg.Nh-1)-(cfg.Nh-1)/2)*cfg.dh_lambda*lambda;
z = ((0:cfg.Nv-1)-(cfg.Nv-1)/2)*cfg.dv_lambda*lambda;
[X,Z] = ndgrid(x,z);
model.family = 'subarray';
model.lambda_m = lambda;
model.Nh = cfg.Nh;
model.Nv = cfg.Nv;
model.flat_positions_m = [X(:),zeros(numel(X),1),Z(:)];
model.x_bounds_m = [min(x) max(x)];
model.z_bounds_m = [min(z) max(z)];
model.aperture_width_m = max(x)-min(x);
model.aperture_height_m = max(z)-min(z);
model.subarray_m = s.m;
model.subarray_n = s.n;
model.subarray_count = s.m*s.n;
model.subarray_Nh = cfg.Nh/s.m;
model.subarray_Nv = cfg.Nv/s.n;
model.subarray_ids = zeros(cfg.Nh*cfg.Nv,1);
model.element_groups = cell(model.subarray_count,1);
model.subarray_centers_m = zeros(model.subarray_count,3);
model.subarray_indices = zeros(model.subarray_count,2);
model.parameter_count = 3*model.subarray_count;
model.parameter_names = cell(1,model.parameter_count);
model.parameter_units = cell(1,model.parameter_count);
model.grids = cell(1,model.parameter_count);
indexGrid = reshape(1:cfg.Nh*cfg.Nv,cfg.Nh,cfg.Nv);
names = {'d','yaw','tilt'};
units = {'lambda_c','deg','deg'};
for j = 1:s.n
    for i = 1:s.m
        b = i+(j-1)*s.m;
        ix = (i-1)*model.subarray_Nh+(1:model.subarray_Nh);
        iz = (j-1)*model.subarray_Nv+(1:model.subarray_Nv);
        indices = indexGrid(ix,iz);
        indices = indices(:);
        model.element_groups{b} = indices;
        model.subarray_ids(indices) = b;
        model.subarray_centers_m(b,:) = mean(model.flat_positions_m(indices,:),1);
        model.subarray_indices(b,:) = [i j];
        for d = 1:3
            p = 3*(b-1)+d;
            model.parameter_names{p} = sprintf('%s_%d_%d',names{d},i,j);
            model.parameter_units{p} = units{d};
            model.grids{p} = grids{d};
        end
    end
end
model.units = model.parameter_units;
model.zero = zeros(1,model.parameter_count);
model.parameter_bounds = zeros(model.parameter_count,2);
for p = 1:model.parameter_count
    model.parameter_bounds(p,:) = [min(model.grids{p}) max(model.grids{p})];
end
model.check_positions_m = model.flat_positions_m;
model.check_grid_size = [cfg.Nh cfg.Nv];
model.constraint_note = 'Rigid blocks; center displacement and yaw/tilt bounds; all element-center pair distances. No panel collision or continuous-surface strain/curvature constraint.';
model.continuous_surface = false;
end

function [geom,valid,reason] = sa_buildGeometry(model,parameters,cfg)
parameters = parameters(:).';
assert(isnumeric(parameters) && isreal(parameters) && numel(parameters)==model.parameter_count && all(isfinite(parameters)),'Invalid subarray parameter vector.');
assert(strcmp(model.family,'subarray'),'Expected a subarray model.');
N = size(model.flat_positions_m,1);
B = model.subarray_count;
blockParameters = reshape(parameters,3,B).';
position = model.flat_positions_m;
rotations = repmat(eye(3),1,1,N);
normal = repmat([0 1 0],N,1);
centers = model.subarray_centers_m;
blockRotations = zeros(3,3,B);
valid = true;
reason = '';
for p = 1:model.parameter_count
    bounds = model.parameter_bounds(p,:);
    tolerance = 1e-10*max(1,max(abs(bounds)));
    if parameters(p)<bounds(1)-tolerance || parameters(p)>bounds(2)+tolerance
        valid = false;
        reason = sa_addReason(reason,['parameter bound: ' model.parameter_names{p}]);
    end
end
for b = 1:B
    indices = model.element_groups{b};
    d = blockParameters(b,1)*model.lambda_m;
    yaw = blockParameters(b,2);
    tilt = blockParameters(b,3);
    c = model.subarray_centers_m(b,:);
    centers(b,:) = c+[0 d 0];
    R = [cosd(yaw) sind(yaw) 0;-sind(yaw) cosd(yaw) 0;0 0 1]*[1 0 0;0 cosd(tilt) -sind(tilt);0 sind(tilt) cosd(tilt)];
    blockRotations(:,:,b) = R;
    % 零转角保留原始x/z字节值，避免中心减加导致平面极限的舍入差异。
    if yaw==0 && tilt==0
        position(indices,2) = model.flat_positions_m(indices,2)+d;
    else
        offset = bsxfun(@minus,model.flat_positions_m(indices,:),c);
        position(indices,:) = bsxfun(@plus,offset*R.',centers(b,:));
    end
    normal(indices,:) = repmat(R(:,2).',numel(indices),1);
    for k = 1:numel(indices)
        rotations(:,:,indices(k)) = R;
    end
end
assert(all(isfinite(position(:))),'Subarray positions are nonfinite.');
geom.positions_m = position;
geom.normals = normal;
geom.R = rotations;
geom.yaw_deg = atan2d(normal(:,1),normal(:,2));
geom.tilt_deg = atan2d(normal(:,3),sqrt(sum(normal(:,1:2).^2,2)));
geom.parameters = parameters;
geom.family = model.family;
geom.flat_positions_m = model.flat_positions_m;
geom.subarray_ids = model.subarray_ids;
geom.element_groups = model.element_groups;
geom.subarray_centers_m = centers;
geom.subarray_flat_centers_m = model.subarray_centers_m;
geom.subarray_parameters = blockParameters;
geom.subarray_R = blockRotations;
geom.subarray_indices = model.subarray_indices;
geom.subarray_m = model.subarray_m;
geom.subarray_n = model.subarray_n;
geom.subarray_Nh = model.subarray_Nh;
geom.subarray_Nv = model.subarray_Nv;
geom.continuous_surface = false;
geom.constraint_note = model.constraint_note;
geom.check_grid_size = model.check_grid_size;
displacement = sqrt(sum((position-model.flat_positions_m).^2,2));
geom.metrics.max_displacement_m = max(displacement);
geom.metrics.max_center_displacement_m = max(abs(blockParameters(:,1)))*model.lambda_m;
geom.metrics.max_normal_angle_deg = max(acosd(max(-1,min(1,normal(:,2)))));
geom.metrics.max_strain = 0; % 每个刚性子阵内部；子阵边界无连续应变定义。
geom.metrics.max_curvature_inv_m = 0; % 各平面子阵内部，不用于限制子阵间转角。
geom.metrics.min_spacing_m = sa_minimumElementSpacing(position);
geom.metrics.constraint_sample_count = N;
geom.metrics.rigid_blocks = true;
if geom.metrics.min_spacing_m<cfg.subarray.min_spacing_lambda*model.lambda_m-1e-10*model.lambda_m
    valid = false;
    reason = sa_addReason(reason,'minimum element spacing');
end
if isfield(cfg,'bs') && isfield(cfg,'users')
    deployment = fimScene('prepare',cfg);
else
    deployment.position_m = [0 0 0];
    deployment.R_global_from_array = eye(3);
    deployment.bs.position_m = deployment.position_m;
    deployment.bs.R_global_from_array = deployment.R_global_from_array;
    deployment.metadata.note = 'Legacy local geometry uses identity deployment.';
end
Rbs = deployment.R_global_from_array;
geom.bs_deployment = deployment;
geom.positions_global_m = bsxfun(@plus,position*Rbs.',deployment.position_m);
geom.flat_positions_global_m = bsxfun(@plus,model.flat_positions_m*Rbs.',deployment.position_m);
geom.normals_global = normal*Rbs.';
geom.subarray_centers_global_m = bsxfun(@plus,centers*Rbs.',deployment.position_m);
geom.R_global = zeros(3,3,N);
for k = 1:N
    geom.R_global(:,:,k) = Rbs*rotations(:,:,k);
end
geom.feasible = valid;
geom.feasibility_reason = reason;
end

function sa_checkPositiveInteger(value,name)
assert(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value) && value>=1 && value==fix(value),['Invalid ' name '.']);
end

function grid = sa_withZero(grid)
assert(isnumeric(grid) && isreal(grid) && isvector(grid) && ~isempty(grid) && all(isfinite(grid(:))),'Subarray grids must be nonempty finite real vectors.');
grid = unique([grid(:).',0]);
end

function d = sa_minimumElementSpacing(position)
d = Inf;
for k = 1:size(position,1)-1
    delta = bsxfun(@minus,position(k+1:end,:),position(k,:));
    d = min(d,min(sqrt(sum(delta.^2,2))));
end
end

function reason = sa_addReason(reason,newReason)
if isempty(reason)
    reason = newReason;
else
    reason = [reason '; ' newReason];
end
end

function model = pc_prepareModel(cfg)
assert(isfield(cfg,'pca') && isstruct(cfg.pca) && isscalar(cfg.pca),'Missing cfg.pca.');
p = cfg.pca;
pc_checkPositiveInteger(cfg.Nh,'Nh');
pc_checkPositiveInteger(cfg.Nv,'Nv');
assert(isnumeric(p.min_spacing_lambda) && isreal(p.min_spacing_lambda) && isscalar(p.min_spacing_lambda) && isfinite(p.min_spacing_lambda) && p.min_spacing_lambda>=0,'Invalid pca.min_spacing_lambda.');
assert(ischar(p.region_mode) && any(strcmp(p.region_mode,{'per_element','shared'})),...
    'pca.region_mode must be per_element or shared.');
grids = {pc_withZero(p.x_grid_lambda),pc_withZero(p.z_grid_lambda)};
lambda = cfg.c0/cfg.fc_Hz;
x = ((0:cfg.Nh-1)-(cfg.Nh-1)/2)*cfg.dh_lambda*lambda;
z = ((0:cfg.Nv-1)-(cfg.Nv-1)/2)*cfg.dv_lambda*lambda;
[X,Z] = ndgrid(x,z);
N = numel(X);
model.family = 'pca';
model.lambda_m = lambda;
model.Nh = cfg.Nh;
model.Nv = cfg.Nv;
model.flat_positions_m = [X(:),zeros(N,1),Z(:)];
model.x_bounds_m = [min(x) max(x)];
model.z_bounds_m = [min(z) max(z)];
model.aperture_width_m = max(x)-min(x);
model.aperture_height_m = max(z)-min(z);
model.region_mode = p.region_mode;
model.shared_x_bounds_lambda = [];
model.shared_z_bounds_lambda = [];
model.parameter_count = 2*N;
model.parameter_names = cell(1,2*N);
model.parameter_units = repmat({'lambda_c'},1,2*N);
model.grids = cell(1,2*N);
model.parameter_bounds = zeros(2*N,2);
model.element_indices = zeros(N,2);
for v = 1:cfg.Nv
    for h = 1:cfg.Nh
        k = h+(v-1)*cfg.Nh;
        model.element_indices(k,:) = [h v];
        model.parameter_names{2*k-1} = sprintf('dx_%d_%d',h,v);
        model.parameter_names{2*k} = sprintf('dz_%d_%d',h,v);
        model.grids{2*k-1} = grids{1};
        model.grids{2*k} = grids{2};
        model.parameter_bounds(2*k-1,:) = [min(grids{1}) max(grids{1})];
        model.parameter_bounds(2*k,:) = [min(grids{2}) max(grids{2})];
    end
end
if strcmp(p.region_mode,'shared')
    bx = pc_resolveBounds(p.x_bounds_lambda,model.x_bounds_m/lambda,'pca.x_bounds_lambda');
    bz = pc_resolveBounds(p.z_bounds_lambda,model.z_bounds_m/lambda,'pca.z_bounds_lambda');
    model.shared_x_bounds_lambda = bx;
    model.shared_z_bounds_lambda = bz;
    % 比较必须允许原始阵列作为零状态；不把它静默投影到新的区域。
    tolerance = 1e-10*max([1,abs(bx),abs(bz)]);
    assert(min(X(:)/lambda)>=bx(1)-tolerance && max(X(:)/lambda)<=bx(2)+tolerance && ...
        min(Z(:)/lambda)>=bz(1)-tolerance && max(Z(:)/lambda)<=bz(2)+tolerance,...
        'Shared PCA region must contain every initial element center so the fixed baseline remains feasible.');
    for k = 1:N
        xb = bx-X(k)/lambda;
        zb = bz-Z(k)/lambda;
        model.parameter_bounds(2*k-1,:) = pc_intersectBounds(model.parameter_bounds(2*k-1,:),xb);
        model.parameter_bounds(2*k,:) = pc_intersectBounds(model.parameter_bounds(2*k,:),zb);
    end
end
% 每行依次为 dx_min,dx_max,dz_min,dz_max；shared 模式已与公共边界求交。
model.displacement_bounds_lambda = [model.parameter_bounds(1:2:end,:),model.parameter_bounds(2:2:end,:)];
model.movement_x_bounds_m = bsxfun(@plus,model.displacement_bounds_lambda(:,1:2)*lambda,X(:));
model.movement_z_bounds_m = bsxfun(@plus,model.displacement_bounds_lambda(:,3:4)*lambda,Z(:));
model.units = model.parameter_units;
model.zero = zeros(1,model.parameter_count);
model.check_positions_m = model.flat_positions_m;
model.check_grid_size = [cfg.Nh cfg.Nv];
model.continuous_surface = false;
model.min_spacing_lambda = p.min_spacing_lambda;
model.constraint_note = ['Independent phase centers translate in the original array-local x-z plane; fixed element frames. ' ...
    'Bounded displacement rectangles, optional shared absolute-position rectangle, and all element-center pair distances. ' ...
    'No mutual coupling, package collision, trajectory collision, or material strain/curvature model.'];
end

function [geom,valid,reason] = pc_buildGeometry(model,parameters,cfg)
assert(strcmp(model.family,'pca'),'Expected a PCA model.');
assert(isnumeric(parameters) && isreal(parameters) && isvector(parameters) && ...
    numel(parameters)==model.parameter_count && all(isfinite(parameters(:))),'Invalid PCA parameter vector.');
parameters = parameters(:).';
N = size(model.flat_positions_m,1);
displacement = reshape(parameters,2,N).';
position = model.flat_positions_m;
position(:,1) = position(:,1)+displacement(:,1)*model.lambda_m;
position(:,3) = position(:,3)+displacement(:,2)*model.lambda_m;
assert(all(isfinite(position(:))),'PCA positions are nonfinite.');
valid = true;
reason = '';
withinBounds = true;
for d = 1:model.parameter_count
    bounds = model.parameter_bounds(d,:);
    tolerance = 1e-10*max(1,max(abs(bounds)));
    if parameters(d)<bounds(1)-tolerance || parameters(d)>bounds(2)+tolerance
        withinBounds = false;
        reason = pc_addReason(reason,['movement bound: ' model.parameter_names{d}]);
    end
end
minimumSpacing = pc_minimumElementSpacing(position);
spacingValid = minimumSpacing>=cfg.pca.min_spacing_lambda*model.lambda_m-1e-10*model.lambda_m;
if ~spacingValid
    reason = pc_addReason(reason,'minimum element spacing');
end
valid = withinBounds && spacingValid;
geom.positions_m = position;
geom.normals = repmat([0 1 0],N,1);
geom.R = repmat(eye(3),1,1,N);
geom.yaw_deg = zeros(N,1);
geom.tilt_deg = zeros(N,1);
geom.parameters = parameters;
geom.pca_displacements_lambda = displacement;
geom.pca_displacements_m = displacement*model.lambda_m;
geom.family = 'pca';
geom.flat_positions_m = model.flat_positions_m;
geom.region_mode = model.region_mode;
geom.shared_x_bounds_lambda = model.shared_x_bounds_lambda;
geom.shared_z_bounds_lambda = model.shared_z_bounds_lambda;
geom.displacement_bounds_lambda = model.displacement_bounds_lambda;
geom.movement_x_bounds_m = model.movement_x_bounds_m;
geom.movement_z_bounds_m = model.movement_z_bounds_m;
geom.continuous_surface = false;
geom.constraint_note = model.constraint_note;
geom.check_grid_size = model.check_grid_size;
geom.metrics.max_displacement_m = max(sqrt(sum((position-model.flat_positions_m).^2,2)));
geom.metrics.max_normal_angle_deg = 0;
geom.metrics.max_strain = 0; % 仅兼容公共报告字段；独立阵元间不定义材料应变。
geom.metrics.max_curvature_inv_m = 0; % PCA 不对任何连续曲面计算曲率。
geom.metrics.min_spacing_m = minimumSpacing;
geom.metrics.constraint_sample_count = N;
geom.metrics.movement_bounds_satisfied = withinBounds;
geom.metrics.minimum_spacing_satisfied = spacingValid;
geom.metrics.original_plane_satisfied = all(position(:,2)==0);
if isfield(cfg,'bs') && isfield(cfg,'users')
    deployment = fimScene('prepare',cfg);
else
    deployment.position_m = [0 0 0];
    deployment.R_global_from_array = eye(3);
    deployment.bs.position_m = deployment.position_m;
    deployment.bs.R_global_from_array = deployment.R_global_from_array;
    deployment.metadata.note = 'Legacy local geometry uses identity deployment.';
end
Rbs = deployment.R_global_from_array;
geom.bs_deployment = deployment;
geom.positions_global_m = bsxfun(@plus,position*Rbs.',deployment.position_m);
geom.flat_positions_global_m = bsxfun(@plus,model.flat_positions_m*Rbs.',deployment.position_m);
geom.normals_global = geom.normals*Rbs.';
geom.R_global = repmat(Rbs,1,1,N);
geom.feasible = valid;
geom.feasibility_reason = reason;
end

function grid = pc_withZero(grid)
assert(isnumeric(grid) && isreal(grid) && isvector(grid) && ~isempty(grid) && all(isfinite(grid(:))),...
    'PCA movement grids must be nonempty finite real vectors.');
grid = unique([grid(:).',0]);
end

function bounds = pc_resolveBounds(bounds,defaultBounds,name)
if isempty(bounds)
    bounds = defaultBounds;
end
assert(isnumeric(bounds) && isreal(bounds) && numel(bounds)==2 && all(isfinite(bounds(:))),['Invalid ' name '.']);
bounds = bounds(:).';
assert(bounds(1)<=bounds(2),[name ' must be ordered [min max].']);
end

function bounds = pc_intersectBounds(a,b)
bounds = [max(a(1),b(1)),min(a(2),b(2))];
end

function pc_checkPositiveInteger(value,name)
assert(isnumeric(value) && isreal(value) && isscalar(value) && isfinite(value) && value>=1 && value==fix(value),['Invalid ' name '.']);
end

function d = pc_minimumElementSpacing(position)
d = Inf;
for k = 1:size(position,1)-1
    delta = bsxfun(@minus,position(k+1:end,:),position(k,:));
    d = min(d,min(sqrt(sum(delta.^2,2))));
end
end

function reason = pc_addReason(reason,newReason)
if isempty(reason)
    reason = newReason;
else
    reason = [reason '; ' newReason];
end
end


function geom = buildOrientationGeometry(flat,states,selection,cfg)
geom=flat; geom.family='orientation'; selection=double(selection(:));
assert(numel(selection)==cfg.N && all(selection>=1 & selection<=states.count),'Invalid pattern-rotation selection.');
angles=[states.yaw(selection),states.tilt(selection)];
geom.parameters=reshape(angles.',1,[]); geom.orientation_parameters_deg=angles;
for n=1:cfg.N
    yaw=angles(n,1); tilt=angles(n,2);
    R=[cosd(yaw) sind(yaw) 0;-sind(yaw) cosd(yaw) 0;0 0 1]*[1 0 0;0 cosd(tilt) -sind(tilt);0 sind(tilt) cosd(tilt)];
    geom.R(:,:,n)=R; geom.normals(n,:)=R(:,2).';
    geom.R_global(:,:,n)=cfg.bs.R_global_from_array*R;
end
geom.normals_global=geom.normals*cfg.bs.R_global_from_array.';
geom.yaw_deg=angles(:,1); geom.tilt_deg=angles(:,2);
geom.metrics.max_displacement_m=0;
geom.metrics.max_normal_angle_deg=max(acosd(max(-1,min(1,geom.normals(:,2)))));
geom.constraint_note='Fixed element positions; independently rotated scalar power patterns; no displacement.';
end
