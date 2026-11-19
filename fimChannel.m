function varargout = fimChannel(action,varargin)
%FIMCHANNEL 信道、阵列和方向图公共模块；各算法作为本文件局部函数。
% prepare(cfg) -> pattern, profile, rawCDL
% realization(cfg,profile) -> realization
% time(realization,time_s) -> 同一子径实现的新时间快照
% series(realization,pattern,geom,cfg,time_s) -> H (N x K x F x T), info
% time_s 是相对于子径初相的时间向量，单位为秒。
% candidates(realization,pattern,states,cfg) -> Hcand (N x K x C x F)
% F 是实际求值频点数；每个阵元的状态 c 在所有 F 个频点上共用。
% geometry(realization,pattern,geom,cfg) -> H (N x K x F), info
% geom.positions_m 和 geom.R 分别描述天面本地位置及阵元本体到天面本地的完整旋转。
% CP 不足的形变候选返回 H=[]、info.cp_feasible=false，供优化器跳过。
% pattern_check(pattern,cfg) -> rotated spherical integrals
% selftest(cfg) -> numerical fixture check labels
switch action
    case 'prepare'
        cfg = varargin{1};
        varargout{1} = loadAndNormalizePattern(cfg);
        [varargout{2},varargout{3}] = loadCDLProfile(cfg);
    case 'realization'
        varargout{1} = buildChannelRealization(varargin{:});
    case 'time'
        varargout{1} = channelTimeSnapshot(varargin{:});
    case 'series'
        [varargout{1},varargout{2}] = buildGeometrySeries(varargin{:});
    case 'candidates'
        varargout{1} = buildFIMCandidateChannels(varargin{:});
    case 'geometry'
        [varargout{1},varargout{2}] = buildGeometryChannel(varargin{:});
    case 'pattern_check'
        varargout{1} = checkPatternRotations(varargin{:});
    otherwise
        error('Unknown fimChannel action: %s',action);
end
end

function pattern = loadAndNormalizePattern(cfg)
%LOADANDNORMALIZEPATTERN 先读完整实测网格，再在功率域完成周期延拓与归一化。
assert(exist(cfg.pattern.file,'file')==2,'Missing antenna pattern TXT: %s',cfg.pattern.file);
phi = 0:359;
theta = 0:180;
[Vgain,headers] = readConfiguredPattern(cfg.pattern.file,phi,theta);
assert(isequal(size(Vgain),[181 360]),'getAntennaPattern must return Theta x Phi grid.');
assert(~any(isnan(Vgain(:))) && ~any(Vgain(:)==Inf),'Invalid measured gain samples.');
G = 10.^(Vgain/10);
[P,T] = meshgrid(phi,theta);
native = [sind(T(:)).*cosd(P(:)),sind(T(:)).*sind(P(:)),cosd(T(:))].';
configuredRotation = cfg.pattern.native_to_body;
actualRotation = configuredRotation;
[~,peakIndex] = max(G(:));
peakNative = native(:,peakIndex);
if strcmp(cfg.pattern.alignment,'peak')
    % 先应用用户的原生坐标标定，再校正主峰方位/仰角到本体 +y。
    % 用当前本体 +z 保持竖直参考；此步骤不代替用户对 native 坐标/roll 的标定。
    % directionFrame 的第三轴是全局 z 在主峰切平面上的投影，极点时用 x。
    preliminaryPeak = configuredRotation*peakNative;
    correction = directionFrame([0;1;0])*directionFrame(preliminaryPeak)';
    actualRotation = correction*configuredRotation;
end
body = actualRotation*native;
if cfg.front_only
    G(reshape(body(2,:)<-1e-12,size(G))) = 0;
end
phi = 0:360;
G = [G,G(:,1)];
rawIntegral = trapz(theta*pi/180,trapz(phi*pi/180,G,2).*sind(theta(:)));
assert(isfinite(rawIntegral) && rawIntegral>0,'Pattern spherical integral must be positive.');
G = G*(4*pi/rawIntegral);
pattern = struct('phi',phi,'theta',theta,'G',G,'native_to_body',actualRotation,'front_only',cfg.front_only);
pattern.configured_native_to_body = configuredRotation;
pattern.alignment = cfg.pattern.alignment;
pattern.source_file = cfg.pattern.file;
pattern.file_headers = headers;
pattern.is_test_pattern = ~isempty(strfind(upper([headers{1} ' ' headers{2}]),'TEST'));
pattern.rawIntegral = rawIntegral;
pattern.scale = 4*pi/rawIntegral;
pattern.integral = trapz(theta*pi/180,trapz(phi*pi/180,G,2).*sind(theta(:)));
assert(abs(pattern.integral-4*pi)<cfg.check.norm_tol*4*pi,'Pattern normalization failed.');
% 预期半功率宽度仅用于对照，不改变读取的方向图或强制拟合。
pattern.diagnostics = patternDiagnostics(pattern);
end

function [gain,headers] = readConfiguredPattern(filename,phi,theta)
% 默认实测文件仍调用原 getAntennaPattern 接口；显式替换路径按相同格式读取。
% 合成 TEST 文件只能通过 cfg.pattern.file 显式选用，不覆盖或回退实测文件。
fid = fopen(filename,'r');
assert(fid~=-1,'Cannot open configured pattern TXT: %s',filename);
cleanup = onCleanup(@() fclose(fid));
headers = {fgetl(fid),fgetl(fid)};
assert(ischar(headers{1}) && ischar(headers{2}),'Pattern TXT must contain two header lines.');
readerDir = fileparts(which('getAntennaPattern'));
expected = fullfile(readerDir,'HW_Airview_Fre6700M_E06_1to6_12.7dB.txt');
if strcmp(filename,expected)
    [gain,~] = getAntennaPattern(phi,theta);
else
    data = fscanf(fid,'%f %f %f %f',[4 inf]).';
    assert(size(data,1)==360*181 && size(data,2)==4,'Expected 360*181 Phi/Theta/gain/phase rows after two headers.');
    [P,T] = ndgrid(0:359,0:180);
    assert(all(abs(data(:,1)-P(:))<1e-8) && all(abs(data(:,2)-T(:))<1e-8),'Pattern order must be Phi-fast and Theta-slow.');
    gain = reshape(data(:,3),360,181).';
end
end

function diagnostics = patternDiagnostics(pattern)
% 主峰切面的最近两侧 -3 dB 交点；窄波束精度仍受原 TXT 采样分辨率限制。
% front_only 和 alignment 后的实际功率图用于诊断，不假定它必有给定宽度。
withoutDuplicate = pattern.G(:,1:end-1);
[peak,idx] = max(withoutDuplicate(:));
[it,ip] = ind2sub(size(withoutDuplicate),idx);
phi = pattern.phi(ip);
theta = pattern.theta(it);
nativePeak = [sind(theta)*cosd(phi);sind(theta)*sind(phi);cosd(theta)];
peakBody = pattern.native_to_body*nativePeak;
az = atan2d(peakBody(1),peakBody(2));
el = atan2d(peakBody(3),sqrt(sum(peakBody(1:2).^2)));
diagnostics.peak_native_phi_theta_deg = [phi theta];
diagnostics.peak_array_az_el_deg = [az el];
diagnostics.peak_gain_dBi = 10*log10(peak);
diagnostics.width_level_dB = -3;
diagnostics.note = 'Nearest -3 dB crossings of principal cuts through the loaded, aligned, front-masked sampled peak; no width fitting.';
offset = 0:0.1:180;
leftAz = az-offset;
rightAz = az+offset;
leftU = [cosd(el)*sind(leftAz);cosd(el)*cosd(leftAz);sind(el)*ones(size(offset))];
rightU = [cosd(el)*sind(rightAz);cosd(el)*cosd(rightAz);sind(el)*ones(size(offset))];
threshold = peak*10^(-3/10);
left = crossingOffset(offset,patternPowerOrientation(leftU,eye(3),pattern),threshold);
right = crossingOffset(offset,patternPowerOrientation(rightU,eye(3),pattern),threshold);
diagnostics.hpbw_az_deg = left+right;
diagnostics.hpbw_az_valid = isfinite(left) && isfinite(right);
leftOffset = offset(offset<=el+90+1e-10);
rightOffset = offset(offset<=90-el+1e-10);
leftEl = el-leftOffset;
rightEl = el+rightOffset;
leftU = [cosd(leftEl)*sind(az);cosd(leftEl)*cosd(az);sind(leftEl)];
rightU = [cosd(rightEl)*sind(az);cosd(rightEl)*cosd(az);sind(rightEl)];
left = crossingOffset(leftOffset,patternPowerOrientation(leftU,eye(3),pattern),threshold);
right = crossingOffset(rightOffset,patternPowerOrientation(rightU,eye(3),pattern),threshold);
diagnostics.hpbw_el_deg = left+right;
diagnostics.hpbw_el_valid = isfinite(left) && isfinite(right);
end

function crossing = crossingOffset(offset,power,threshold)
crossing = NaN;
i = find(power<=threshold,1,'first');
if isempty(i) || i==1
    return
end
fraction = (power(i-1)-threshold)/(power(i-1)-power(i));
crossing = offset(i-1)+fraction*(offset(i)-offset(i-1));
end

function [profile,raw] = loadCDLProfile(cfg)
%LOADCDLPROFILE 优先读取已有 MAT，否则调用配置的 CDL 簇模板生成器。
if ~isempty(cfg.cdl.profile_mat)
    raw = readCDLProfileFile(cfg.cdl.profile_mat);
else
    generator = cfg.cdl.generator;
    assert(~isempty(which(generator)),'CDL generator function not found: %s',generator);
    raw = feval(generator,['NR-' cfg.cdl.profile],'CenterFrequency',cfg.fc_Hz,cfg.cdl.extra_args{:});
end
required = {'pow','taus'};
if strcmp(cfg.cdl.bs_angle_side,'arrival')
    af = 'AoA';
    ef = 'EoA';
else
    af = 'AoD';
    ef = 'EoD';
end
required = [required,{af,ef}];
for i = 1:numel(required)
    assert(isfield(raw,required{i}),['Missing CDL field: ' required{i}]);
    v = raw.(required{i});
    assert(isnumeric(v) && isreal(v) && isvector(v) && all(isfinite(v(:))),['Invalid CDL field: ' required{i}]);
end
p = double(raw.pow(:));
tau = double(raw.taus(:));
az = raw.(af);
ea = raw.(ef);
az = double(az(:));
ea = double(ea(:));
assert(numel(p)==numel(tau) && numel(p)==numel(az) && numel(p)==numel(ea),'CDL field lengths differ.');
assert(all(p>=0) && sum(p)>0 && all(tau>=0),'pow must be nonnegative linear powers, taus must be seconds.');
assert(all(abs(az)<=2*pi+1e-8),'AoA/AoD must be radians, not degrees.');
if strcmp(cfg.cdl.eoa_convention,'zenith')
    assert(all(ea>=-1e-8 & ea<=pi+1e-8),'Zenith angles must lie in [0,pi].');
    el = pi/2-ea;
else
    assert(all(abs(ea)<=pi/2+1e-8),'EoA/EoD is outside elevation range; check cdl.eoa_convention and the input angles.');
    el = ea;
end
profile.p = p/sum(p);
profile.tau = tau;
profile.U = [cos(el).*cos(az),cos(el).*sin(az),sin(el)].';
profile.numPaths = numel(p);
modelName = cfg.cdl.profile;
if isfield(raw,'name')
    modelName = char(raw.name);
end
profile.model = [modelName ' scalar cluster approximation'];
profile.bs_angle_side = cfg.cdl.bs_angle_side;
profile.eoa_convention = cfg.cdl.eoa_convention;
if isfield(raw,'numClusters')
    assert(raw.numClusters==numel(p),'Expected one angle/power entry per cluster.');
end
if strcmp(cfg.cdl.model,'rays')
    % Validate the ray template at input preparation without advancing the RNG.
    options = cfg.cdl;
    options.ray_coupling = 'aligned';
    expandCDLRays(raw,options);
    profile.ray_template = raw;
    profile.model = [modelName ' 3GPP ray-level CDL'];
end
end

function raw = readCDLProfileFile(profileFile)
% 局部作用域内无返回值 load，兼容不支持 s=load(...) 的 MLang。
% 外部 MAT 的变量不会覆盖主信道函数的 cfg/profile 等已有状态。
load(profileFile);
if exist('rawCDL','var')==1
    raw = rawCDL;
elseif exist('result','var')==1 && isstruct(result) && isfield(result,'rawCDL')
    raw = result.rawCDL;
else
    error('MAT must contain rawCDL or result.rawCDL.');
end
end

function [U,info] = rotateCDLAnglesToUser(U0,p,target,cfg)
%ROTATECDLANGLESTOUSER 单个三维刚性旋转保持所有簇之间的夹角。
p = p(:);
target = target(:)/norm(target);
assert(size(U0,1)==3 && size(U0,2)==numel(p),'Profile shape mismatch.');
switch cfg.cdl.center_mode
    case 'power_mean'
        reference = U0*p;
        assert(norm(reference)>1e-8,'CDL power-weighted direction is ambiguous; configure explicit reference.');
        reference = reference/norm(reference);
    case 'explicit'
        a = cfg.cdl.reference_az_el_deg;
        assert(numel(a)==2 && all(isfinite(a)) && abs(a(2))<=90,'Set reference [standard az, elevation].');
        reference = [cosd(a(2))*cosd(a(1));cosd(a(2))*sind(a(1));sind(a(2))];
    otherwise
        error('Unknown CDL center_mode.');
end
% 先用一个正交局部框架移走 profile 中心，再放到目标；不用 az/el 直接相加。
R = directionFrame(target)*directionFrame(reference)';
U = R*U0;
info.R = R;
info.reference = reference;
info.target = target;
end

function F = directionFrame(u)
% 右手框架的第一列指向中心，使用当前参考坐标的 +z 定义其余轴，极点改用 +x。
up = [0;0;1];
if abs(up'*u)>1-1e-10
    up = [1;0;0];
end
side = cross(up,u);
side = side/norm(side);
F = [u,side,cross(u,side)];
end

function R = elementRotation(yaw,tilt)
%ELEMENTROTATION 阵元本体到 BS 阵列局部坐标的主动旋转：Rz(-yaw)*Rx(tilt)。
% R*[0;1;0] = [cos(tilt)*sin(yaw); cos(tilt)*cos(yaw); sin(tilt)].
c = cosd(yaw);
s = sind(yaw);
ct = cosd(tilt);
st = sind(tilt);
R = [c s 0; -s c 0; 0 0 1]*[1 0 0; 0 ct -st; 0 st ct];
end

function G = patternPowerRotated(U,yaw,tilt,pattern)
%PATTERNPOWERROTATED 3D 逆旋转后做 pairwise 插值；U 为 3 x L。
G = patternPowerOrientation(U,elementRotation(yaw,tilt),pattern);
end

function G = patternPowerOrientation(U,R,pattern)
%PATTERNPOWERORIENTATION 以完整局部坐标系旋转实测方向图，保留 roll。
% U 已在 BS 阵列坐标：先到阵元本体，再到 TXT 原生球坐标，不能再减 BS 下倾。
local = R'*U;
native = pattern.native_to_body'*local;
phi = mod(atan2d(native(2,:),native(1,:)),360);
theta = acosd(max(-1,min(1,native(3,:))));
G = interp2(pattern.phi,pattern.theta,pattern.G,phi,theta,'linear');
assert(all(isfinite(G)),'Pattern query returned invalid gain.');
G = max(G,0);
if pattern.front_only
    G(local(2,:)<-1e-12) = 0;
end
% front_only 的半球随阵元方向图一起旋转；不在全局坐标二次裁切。
end

function [A,pos_lambda] = arrayResponseUPA(U,cfg,frequencyOffset_Hz)
%ARRAYRESPONSEUPA U: 3 x L，方向从阵列指向来波源；A: N x L。
% pos_lambda 始终以载频波长表示；宽带时只改变空间相位，不改变实际位置。
% 不除 sqrt(N)，保留物理接收阵列增益。
if nargin<3
    frequencyOffset_Hz = 0;
end
assert(size(U,1)==3 && all(isfinite(U(:))),'U must be 3 x L and finite.');
assert(max(abs(sum(U.^2,1)-1))<1e-8,'U columns must be unit vectors.');
[X,Z] = ndgrid(((0:cfg.Nh-1)-(cfg.Nh-1)/2)*cfg.dh_lambda,((0:cfg.Nv-1)-(cfg.Nv-1)/2)*cfg.dv_lambda);
pos_lambda = [X(:),zeros(numel(X),1),Z(:)];
ratio = 1;
if cfg.ofdm.frequency_dependent_steering
    ratio = (cfg.fc_Hz+frequencyOffset_Hz)/cfg.fc_Hz;
end
A = exp(1i*2*pi*ratio*(pos_lambda*U));
end

function realization = buildChannelRealization(cfg,profile)
%BUILDCHANNELREALIZATION 一次抽样，fixed/FIM、所有频点和 SNR 共用。
% 每条路径的时延和系数跨频点共用。子径模式另保存两分量入射场。
% 部署位置先经 BS 的全局姿态变换到阵列局部坐标，机械形变仍在局部坐标。
% 假定完美上行定时对齐：每 UE 的公共距离 d/c 已补偿，不加入 CDL 超额时延/CP。
% BS 的全局平移对应的参考载频相位吸收入 alpha；不得重复加入局部空间相位。
scene = fimScene('users',cfg);
realization.userScene = scene;
realization.bs_frame = scene.bs;
realization.coordinate_frame = 'BS array local: +x horizontal, +y front, +z up';
realization.timing_note = 'Perfect uplink timing alignment removes each UE common d/c; tau contains CDL excess delay only.';
realization.users.az = scene.azimuth_array_deg(:);
realization.users.el = scene.elevation_array_deg(:);
realization.users.U = scene.direction_array;
assert(isequal(size(realization.users.U),[3 cfg.K]),'Scene must provide local BS-frame 3 x K directions.');
assert(numel(scene.pathloss_power_gain)==cfg.K && all(isfinite(scene.pathloss_power_gain(:))) && all(scene.pathloss_power_gain(:)>=0),'Invalid scene path-loss power gains.');
realization.paths = cell(cfg.K,1);
rayModel = isfield(profile,'ray_template');
if rayModel
    realization.cdl_time_s = cfg.cdl.time_s;
end
Rbs = scene.bs.R_global_from_array;
earliest = Inf;
latest = -Inf;
for k = 1:cfg.K
    if rayModel
        path = buildCDLRayPath(profile.ray_template,scene,k,cfg);
    else
        if ~cfg.cdl.align_to_ue
            uGlobal = profile.U;
            u = Rbs'*uGlobal;
            rotation = struct('reference_frame','global_template','R_global',eye(3), ...
                'R_array',Rbs','R',Rbs','reference',[],'target',[]);
        elseif strcmp(cfg.users.mode,'angles')
            % 旧角度抽样模式以阵列局部坐标定义簇内竖直参考。
            [u,rotation] = rotateCDLAnglesToUser(profile.U,profile.p,realization.users.U(:,k),cfg);
            uGlobal = Rbs*u;
            rotation.reference_frame = 'array_local_legacy_angles';
            rotation.R_array = rotation.R;
            rotation.R_global = Rbs*rotation.R;
        else
            [uGlobal,rotation] = rotateCDLAnglesToUser(profile.U,profile.p,scene.direction_global(:,k),cfg);
            u = Rbs'*uGlobal;
            rotation.reference_frame = 'global_environment';
            rotation.R_global = rotation.R;
            rotation.R_array = Rbs'*rotation.R;
            rotation.R = rotation.R_array;
        end
        rotation.target_global = scene.direction_global(:,k);
        rotation.target_array = scene.direction_array(:,k);
        rotation.target = rotation.target_array;
        alpha = sqrt(profile.p(:)/2).*(randn(profile.numPaths,1)+1i*randn(profile.numPaths,1));
        alpha = sqrt(scene.pathloss_power_gain(k))*alpha;
        path = struct('U',u,'U_local',u,'U_global',uGlobal,'alpha',alpha, ...
            'tau',profile.tau(:),'rotation',rotation,'pathloss_power_gain',scene.pathloss_power_gain(k));
    end
    realization.paths{k} = path;
    [~,pos_lambda] = arrayResponseUPA(path.U,cfg);
    % h_n(f) 中的总时延为 tau_l - r_n'*u_l/c0。
    % 因 r_n=pos_lambda*c0/fc，几何时延项可写为 pos_lambda*u/fc。
    % 检查所有 UE、阵元、非零功率路径的共同 FFT 窗口所需时延跨度。
    % 即使关闭频变 steering，仍保留物理阵列几何时延的保守 CP 约束。
    activePaths = abs(path.alpha(:))>0;
    assert(any(activePaths),'At least one positive-power CDL path is required.');
    effectiveDelay = bsxfun(@minus,path.tau(activePaths).',pos_lambda*path.U(:,activePaths)/cfg.fc_Hz);
    earliest = min(earliest,min(effectiveDelay(:)));
    latest = max(latest,max(effectiveDelay(:)));
end
realization.delay_min_s = earliest;
realization.delay_max_s = latest;
realization.delay_span_s = latest-earliest;
realization.array_positions_m = pos_lambda*(cfg.c0/cfg.fc_Hz);
if cfg.ofdm.enabled && realization.delay_span_s>cfg.ofdm.cp_duration_s+1e-15
    error('fimChannel:CPTooShort','CP %.6g us is shorter than common FFT-window delay span %.6g us. Increase cfg.ofdm.cp_ratio (ratio mode) or cp_samples (samples mode), or reduce channel delay spread; this model does not simulate ISI/ICI.',cfg.ofdm.cp_duration_s*1e6,realization.delay_span_s*1e6);
end
% 不按 realization 或频点归一化 alpha/H；保留频率选择性深衰落。
end

function path = buildCDLRayPath(raw,scene,k,cfg)
% The receive and transmit polarization bases both use global coordinates.
rays = expandCDLRays(raw,cfg.cdl);
Rbs = scene.bs.R_global_from_array;
rxRotation = eye(3);
txRotation = eye(3);
if cfg.cdl.align_to_ue
    target = scene.direction_global(:,k);
    arrivalReference = rays.referenceArrival;
    if ~any(rays.isLOS) && strcmp(cfg.cdl.center_mode,'explicit')
        a = cfg.cdl.reference_az_el_deg;
        arrivalReference = [cosd(a(2))*cosd(a(1)); cosd(a(2))*sind(a(1)); sind(a(2))];
    end
    rxRotation = directionFrame(target)*directionFrame(arrivalReference)';
    txRotation = directionFrame(-target)*directionFrame(rays.referenceDeparture)';
end
uGlobal = rxRotation*rays.Uarrival;
txGlobal = txRotation*rays.Udeparture;
[incident0,~,doppler_Hz] = cdlRayCoefficients(rays,txGlobal,eye(3),cfg.cdl.tx_slant_deg, ...
    cfg.cdl.velocity_mps,0,cfg.fc_Hz,cfg.c0,scene.distance_3d_m(k));
pathloss = scene.pathloss_power_gain(k);
u = Rbs'*uGlobal;
path = struct('U',u,'U_local',u,'U_global',uGlobal,'tau',rays.tau, ...
    'alpha',sqrt(pathloss*rays.power),'pathloss_power_gain',pathloss, ...
    'incident0',sqrt(pathloss)*incident0,'doppler_Hz',doppler_Hz(:).', ...
    'R_global_from_array',Rbs,'rx_slant_deg',cfg.cdl.rx_slant_deg, ...
    'Udeparture_global',txGlobal,'isLOS',rays.isLOS,'cluster',rays.cluster);
path.incidentField = fieldAtTime(path,cfg.cdl.time_s);
end

function realization = channelTimeSnapshot(realization,time_s)
validateattributes(time_s,{'numeric'},{'real','finite','scalar','nonnegative'});
assert(isfield(realization,'cdl_time_s'),'Time snapshots require a ray-level CDL realization.');
for k = 1:numel(realization.paths)
    path = realization.paths{k};
    path.incidentField = fieldAtTime(path,time_s);
    realization.paths{k} = path;
end
realization.cdl_time_s = time_s;
end

function field = fieldAtTime(path,time_s)
field = path.incident0.*exp(1i*2*pi*time_s*path.doppler_Hz);
end

function a = pathAmplitude(path,R,pattern)
gain = sqrt(patternPowerOrientation(path.U,R,pattern));
if isfield(path,'incidentField')
    polarization = cdlPolarizationField(path.U_global,path.R_global_from_array*R,path.rx_slant_deg);
    a = gain.*sum(polarization.*path.incidentField,1);
else
    a = gain.*path.alpha(:).';
end
end

function [H,info] = buildGeometrySeries(realization,pattern,geom,cfg,time_s)
validateattributes(time_s,{'numeric'},{'real','finite','vector','nonempty','nonnegative'}, ...
    'fimChannel/series','time_s');
assert(isfield(realization,'cdl_time_s'),'Time series require a ray-level CDL realization.');
time_s = reshape(full(double(time_s)),1,[]);
[H,info] = buildGeometryChannel(realization,pattern,geom,cfg,time_s);
info.time_s = time_s;
end

function [H,info] = buildGeometryChannel(realization,pattern,geom,cfg,time_s)
%BUILDGEOMETRYCHANNEL 真实位置和局部完整朝向耦合的宽带信道。
% 本函数不改变 alpha/U/tau，不按形状或频点重新归一化任何信道。
% 方向图功率增益只计算一次，之后在所有频点共用。
positions = geom.positions_m;
rotations = geom.R;
assert(isnumeric(positions) && isreal(positions) && isequal(size(positions),[cfg.N 3]) && all(isfinite(positions(:))),'geom.positions_m must be finite real N x 3.');
assert(isnumeric(rotations) && isreal(rotations) && size(rotations,1)==3 && size(rotations,2)==3 && size(rotations,3)==cfg.N && numel(rotations)==9*cfg.N && all(isfinite(rotations(:))),'geom.R must be finite real 3 x 3 x N.');
for n = 1:cfg.N
    R = rotations(:,:,n);
    assert(norm(R'*R-eye(3),'fro')<1e-9 && abs(det(R)-1)<1e-9,'Each geom.R must be a proper orthogonal rotation.');
end
% 元素全局朝向 Rbs*R_n；实测方向图随 BS 下倾与局部 FIM 旋转一起转动。
% 计算中采用等价的 R_n'*(Rbs'*u_global)，避免再次扣减下倾角。
if isfield(realization,'bs_frame')
    bsFrame = realization.bs_frame;
else
    % 人工/旧版局部信道夹具可无 userScene，部署矩阵只作元数据记录。
    deployment = fimScene('prepare',cfg);
    bsFrame = deployment.bs;
end
Rbs = bsFrame.R_global_from_array;
info.coordinate_frame = 'Channel positions and directions evaluated in BS array-local coordinates';
info.bs_position_m = bsFrame.position_m;
info.R_global_from_array = Rbs;
info.positions_global_m = bsxfun(@plus,positions*Rbs.',bsFrame.position_m);
info.R_global_from_element = zeros(3,3,cfg.N);
info.normals_global = zeros(cfg.N,3);
for n = 1:cfg.N
    info.R_global_from_element(:,:,n) = Rbs*rotations(:,:,n);
    normalGlobal = Rbs*rotations(:,2,n);
    info.normals_global(n,:) = normalGlobal.';
end
assert(iscell(realization.paths) && numel(realization.paths)==cfg.K,'One path structure per UE is required.');
earliest = Inf;
latest = -Inf;
for k = 1:cfg.K
    path = realization.paths{k};
    L = numel(path.alpha);
    assert(L>=1,'Each UE must have at least one path entry.');
    assert(isnumeric(path.U) && isreal(path.U) && isequal(size(path.U),[3 L]) && all(isfinite(path.U(:))),'Path directions must be finite real 3 x L.');
    assert(max(abs(sum(path.U.^2,1)-1))<1e-8,'Path directions must be unit vectors.');
    assert(isnumeric(path.alpha) && isvector(path.alpha) && all(isfinite(path.alpha(:))) && isnumeric(path.tau) && isreal(path.tau) && numel(path.tau)==L && all(isfinite(path.tau(:))),'Invalid path coefficients or reference delays.');
    % 零功率路径 alpha=0 不影响 CP；不依赖姿态遮挡来放松 CP 约束。
    active = abs(path.alpha(:))>0;
    if any(active)
        projection_s = positions*path.U(:,active)/cfg.c0;
        referenceDelay = path.tau(:);
        effectiveDelay = bsxfun(@minus,referenceDelay(active).',projection_s);
        earliest = min(earliest,min(effectiveDelay(:)));
        latest = max(latest,max(effectiveDelay(:)));
    end
end
if isinf(earliest)
    earliest = 0;
    latest = 0;
end
info.delay_min_s = earliest;
info.delay_max_s = latest;
info.delay_span_s = latest-earliest;
% 即使关闭频变 steering，仍采用真实几何时延作保守 CP 检查。
info.cp_feasible = ~cfg.ofdm.enabled || info.delay_span_s<=cfg.ofdm.cp_duration_s+1e-15;
H = [];
if ~info.cp_feasible
    return
end
F = cfg.ofdm.num_tones;
isSeries = nargin>=5;
T = 1;
if isSeries
    T = numel(time_s);
end
H = complex(zeros(cfg.N,cfg.K,F,T));
for k = 1:cfg.K
    path = realization.paths{k};
    if isSeries
        % Project the initial incident field once. Doppler is scalar per ray
        % and can be applied after polarization projection at every time.
        path.incidentField = path.incident0;
        timePhase = exp(1i*2*pi*(path.doppler_Hz(:)*time_s));
    end
    L = numel(path.alpha);
    amplitudes = complex(zeros(cfg.N,L));
    for n = 1:cfg.N
        amplitudes(n,:) = pathAmplitude(path,rotations(:,:,n),pattern);
    end
    projection_s = positions*path.U/cfg.c0;
    for f = 1:F
        nu = cfg.ofdm.freq_offsets_Hz(f);
        spatialFrequency = cfg.fc_Hz;
        if cfg.ofdm.frequency_dependent_steering
            spatialFrequency = cfg.fc_Hz+nu;
        end
        spatialPhase = exp(1i*2*pi*spatialFrequency*projection_s);
        % alpha 已含参考载频相位；这里时延因子仅使用频率偏移 nu。
        % 不能再把几何时延加入 tau，否则位置频率相位会重复计算。
        delayPhase = exp(-1i*2*pi*nu*path.tau(:));
        weightedPaths = spatialPhase.*amplitudes;
        if isSeries
            channel = weightedPaths*(delayPhase.*timePhase);
            H(:,k,f,:) = reshape(channel,cfg.N,1,1,T);
        else
            H(:,k,f) = weightedPaths*delayPhase;
        end
    end
end
assert(all(isfinite(H(:))),'Geometry channel contains nonfinite coefficients.');
end

function values = checkPatternRotations(pattern,cfg)
[P,T] = meshgrid(0:360,0:180);
U = [sind(T(:)).*cosd(P(:)),sind(T(:)).*sind(P(:)),cosd(T(:))].';
angles = [0 0;30 20;-30 -20];
values = zeros(3,1);
for i = 1:3
    G = reshape(patternPowerRotated(U,angles(i,1),angles(i,2),pattern),size(P));
    values(i) = trapz((0:180)*pi/180,trapz((0:360)*pi/180,G,2).*sind((0:180)'));
end
assert(max(abs(values/(4*pi)-1))<0.02,'Rotated integral mismatch >2%%; inspect mapping, grid, and front cut.');
end

function Hcand = buildFIMCandidateChannels(realization,pattern,states,cfg)
%BUILDFIMCANDIDATECHANNELS Hcand(n,k,c,f)：阵元 n、UE k、状态 c、频点 f。
% 不允许每个子载波独立选择 FIM 角度；优化器只对每个阵元选一个 c_n。
F = cfg.ofdm.num_tones;
Hcand = complex(zeros(cfg.N,cfg.K,states.count,F));
for k = 1:cfg.K
    path = realization.paths{k};
    [A0,pos_lambda] = arrayResponseUPA(path.U,cfg);
    phaseProjection = pos_lambda*path.U;
    L = numel(path.alpha);
    % 在所有频点复用每个姿态的路径幅度。
    amplitude = complex(zeros(L,states.count));
    for c = 1:states.count
        R = elementRotation(states.yaw(c),states.tilt(c));
        amplitude(:,c) = pathAmplitude(path,R,pattern).';
    end
    for f = 1:F
        nu = cfg.ofdm.freq_offsets_Hz(f);
        if cfg.ofdm.frequency_dependent_steering
            A = exp(1i*2*pi*((cfg.fc_Hz+nu)/cfg.fc_Hz)*phaseProjection);
        else
            A = A0;
        end
        delayPhase = exp(-1i*2*pi*nu*path.tau(:));
        channels = A*bsxfun(@times,delayPhase,amplitude);
        Hcand(:,k,:,f) = reshape(channels,[cfg.N 1 states.count]);
    end
end
assert(all(isfinite(Hcand(:))),'Invalid candidate channel.');
end

