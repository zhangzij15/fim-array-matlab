function rays = expandCDLRays(raw,options)
%EXPANDCDLRAYS Expand a terrestrial CDL cluster template into 3GPP rays.
% TR 38.901 V18.0.0 §7.7.1 and Table 7.5-3 define 20 rays per diffuse row.
% Directions are 3-by-L; other ray fields are L-by-1. XPR is a linear ratio.
% angle_scaling_deg is [desired means; desired spreads] for AoD/AoA/ZoD/ZoA.
% [NaN; NaN] skips the corresponding angle (angles and spreads in degrees).

if nargin<2
    options = struct();
end
if ~isfield(options,'eoa_convention'), options.eoa_convention = 'elevation'; end
if ~isfield(options,'ray_coupling'), options.ray_coupling = 'random'; end
if ~isfield(options,'k_factor_dB'), options.k_factor_dB = []; end
if ~isfield(options,'angle_scaling_deg'), options.angle_scaling_deg = []; end
convention = validatestring(options.eoa_convention,{'elevation','zenith'});
coupling = validatestring(options.ray_coupling,{'random','aligned'});

assert(isfield(raw,'name'),'FIM:MissingCDLProfileName', ...
    'Ray expansion requires rawCDL.name to identify a terrestrial CDL-A through CDL-E template.');
profile = upper(char(raw.name));
if startsWith(profile,'NR-'), profile = profile(4:end); end
names = {'CDL-A','CDL-B','CDL-C','CDL-D','CDL-E'};
profile = validatestring(profile,names,mfilename,'rawCDL.name');
[tablePaths,tableXPR_dB] = cdlProfileTable(profile);
clusterCount = size(tablePaths,1);
fields = {'pow','taus','AoD','AoA','EoD','EoA'};
values = zeros(clusterCount,numel(fields));
for j = 1:numel(fields)
    field = fields{j};
    assert(isfield(raw,field),'FIM:MissingCDLRayField', ...
        'Ray expansion requires rawCDL.%s.',field);
    validateattributes(raw.(field),{'numeric'}, ...
        {'real','finite','vector','numel',clusterCount},mfilename,['rawCDL.' field]);
    values(:,j) = double(raw.(field)(:));
end
if isfield(raw,'numClusters')
    validateattributes(raw.numClusters,{'numeric'}, ...
        {'real','finite','scalar','integer'},mfilename,'rawCDL.numClusters');
    assert(raw.numClusters==clusterCount,'FIM:CDLClusterCountMismatch', ...
        'rawCDL.numClusters does not match the template row count.');
end
p = values(:,1);
tau = values(:,2);
assert(all(p>=0) && isfinite(sum(p)) && sum(p)>0 && all(tau>=0), ...
    'FIM:InvalidCDLRayPowersDelays','CDL powers must be nonnegative and delays must be nonnegative seconds.');
p = p/sum(p);
clusterAngles = rad2deg(values(:,3:6));
if strcmp(convention,'elevation')
    clusterAngles(:,3:4) = 90-clusterAngles(:,3:4);
end
assert(all(abs(clusterAngles(:,1:2))<=360+1e-8,'all') && ...
    all(clusterAngles(:,3:4)>=-1e-8 & clusterAngles(:,3:4)<=180+1e-8,'all'), ...
    'FIM:InvalidCDLRayAngles','Raw azimuth and elevation or zenith angles must be in radians.');

losProfile = any(strcmp(profile,{'CDL-D','CDL-E'}));
if losProfile
    % The first two table rows are specular and diffuse parts of the same cluster.
    azDifference = mod(clusterAngles(1,1:2)-clusterAngles(2,1:2)+180,360)-180;
    zenDifference = clusterAngles(1,3:4)-clusterAngles(2,3:4);
    assert(p(1)>0 && abs(tau(1)-tau(2))<=1e-15 && ...
        all(abs([azDifference zenDifference])<=1e-8), ...
        'FIM:InvalidCDLSpecular','CDL-D/E must start with a specular row and its co-located diffuse row.');
end
if ~isempty(options.k_factor_dB)
    assert(losProfile,'FIM:CDLKFactorNeedsLOS','K-factor scaling applies to CDL-D and CDL-E.');
    [p,tau] = adjustKFactor(p,tau,raw,options.k_factor_dB);
end

if isfield(raw,'XPR')
    validateattributes(raw.XPR,{'numeric'},{'real','finite','vector','nonempty'}, ...
        mfilename,'rawCDL.XPR');
    assert(isscalar(raw.XPR) || numel(raw.XPR)==clusterCount, ...
        'FIM:InvalidCDLXPR','rawCDL.XPR must be scalar or one value per template row.');
    xpr = 10.^(double(raw.XPR(:))/10);
else
    xpr = 10^(tableXPR_dB/10);
end
assert(all(isfinite(xpr) & xpr>0),'FIM:InvalidCDLXPR', ...
    'CDL XPR must convert to a finite positive linear ratio.');
if isscalar(xpr), xpr = repmat(xpr,clusterCount,1); end

% Normalized offsets (Table 7.5-3) and the four per-cluster spreads
% [cASD cASA cZSD cZSA] from Tables 7.7.1-1 through 7.7.1-5.
offset = [0.0447; 0.1413; 0.2492; 0.3715; 0.5129; ...
    0.6797; 0.8844; 1.1481; 1.5195; 2.1551];
offset = reshape([offset.'; -offset.'],[],1);
spreads = [5 11 3 3; 10 22 3 7; 2 15 3 7; 5 8 3 3; 5 11 3 7];
offset_deg = offset*spreads(strcmp(profile,names),:);
rayCounts = repmat(numel(offset),clusterCount,1);
if losProfile, rayCounts(1) = 1; end
rays.cluster = repelem((1:clusterCount).',rayCounts);
rays.power = p(rays.cluster)./rayCounts(rays.cluster);
rays.tau = tau(rays.cluster);
rays.xpr = xpr(rays.cluster);
rays.isLOS = losProfile & rays.cluster==1;
rayAngles = clusterAngles(rays.cluster,:);
for cluster = (1+losProfile):clusterCount
    first = 1 + (cluster-1)*numel(offset) - losProfile*(numel(offset)-1);
    rows = first:(first+numel(offset)-1);
    increments = offset_deg;
    if strcmp(coupling,'random')
        % Holding AoA fixed and independently permuting the other lists
        % produces the three independent pairings in §7.7.1 Step 2.
        for j = [1 3 4]
            increments(:,j) = increments(randperm(numel(offset)),j);
        end
    end
    rayAngles(rows,:) = rayAngles(rows,:) + increments;
end
if ~isempty(options.angle_scaling_deg)
    [rayAngles,clusterAngles] = scaleAngles(rayAngles,clusterAngles,rays.power,options.angle_scaling_deg);
end
rays.Udeparture = sphericalDirections(rayAngles(:,1),rayAngles(:,3));
rays.Uarrival = sphericalDirections(rayAngles(:,2),rayAngles(:,4));
if losProfile
    referencePower = zeros(size(p));
    referencePower(1) = 1;
else
    referencePower = p;
end
rays.referenceDeparture = directionCenter(clusterAngles(:,1),clusterAngles(:,3),referencePower);
rays.referenceArrival = directionCenter(clusterAngles(:,2),clusterAngles(:,4),referencePower);
end

function [p,tau] = adjustKFactor(p,tau,raw,targetK_dB)
% §7.7.6 scales the sum of diffuse powers relative to the specular path.
validateattributes(targetK_dB,{'numeric'},{'real','finite','scalar'},mfilename,'k_factor_dB');
assert(sum(p(2:end))>0,'FIM:CDLKFactorNeedsDiffuse','K-factor scaling requires diffuse power.');
targetDS = rmsDelay(p,tau);
if isfield(raw,'DS')
    validateattributes(raw.DS,{'numeric'},{'real','finite','scalar','nonnegative'},mfilename,'rawCDL.DS');
    targetDS = double(raw.DS);
end
p(2:end) = p(2:end)./sum(p(2:end));
% Work with the smaller of K and 1/K to avoid overflow.
ratio = 10^(-abs(double(targetK_dB))/10);
assert(ratio>0,'FIM:InvalidCDLKFactor','The requested K-factor is numerically unrepresentable.');
if targetK_dB>=0
    p(1) = 1/(1+ratio);
    p(2:end) = p(2:end)*(ratio/(1+ratio));
else
    p(1) = ratio/(1+ratio);
    p(2:end) = p(2:end)/(1+ratio);
end
if targetDS==0
    tau(:) = 0;
else
    scaledDS = rmsDelay(p,tau);
    assert(scaledDS>0 && isfinite(scaledDS),'FIM:InvalidCDLDelaySpread', ...
        'K-factor scaling requires a positive delay spread for a nonzero target.');
    tau = tau*(targetDS/scaledDS);
end
end

function ds = rmsDelay(p,tau)
ds = sqrt(sum(p.*(tau-sum(p.*tau)).^2));
end

function [angles,centers] = scaleAngles(angles,centers,p,target)
% §7.7.5.1 and Annex A use circular moments of the expanded ray angles.
validateattributes(target,{'numeric'},{'real','size',[2 4]},mfilename,'angle_scaling_deg');
for j = 1:4
    if all(isnan(target(:,j))), continue; end
    validateattributes(target(:,j),{'numeric'},{'finite'},mfilename,'angle_scaling_deg column');
    assert(target(2,j)>=0,'FIM:InvalidCDLAngleScaling','Target angular spreads must be nonnegative.');
    if j>=3
        assert(target(1,j)>=0 && target(1,j)<=180,'FIM:InvalidCDLAngleScaling', ...
            'Target zenith means must be between 0 and 180 degrees.');
    end
    moment = sum(p.*exp(1i*deg2rad(angles(:,j))))/sum(p);
    spread = rad2deg(sqrt(-2*log(min(1,abs(moment)))));
    assert(abs(moment)>eps && spread>0,'FIM:InvalidCDLAngleScaling', ...
        'The input angular mean or spread is undefined for column %d.',j);
    mean_deg = rad2deg(angle(moment));
    if j>=3, mean_deg = mod(mean_deg,360); end
    values = [angles(:,j); centers(:,j)]-mean_deg;
    if j<=2, values = mod(values+180,360)-180; end
    scaled = target(1,j) + (target(2,j)/spread)*values;
    if j<=2
        scaled = mod(scaled,360);
    else
        scaled = min(180,max(0,scaled));
    end
    angles(:,j) = scaled(1:size(angles,1));
    centers(:,j) = scaled(size(angles,1)+1:end);
end
end

function U = sphericalDirections(az_deg,zen_deg)
az = az_deg(:).';
zen = zen_deg(:).';
U = [sind(zen).*cosd(az);sind(zen).*sind(az);cosd(zen)];
end

function reference = directionCenter(az_deg,zen_deg,p)
reference = sphericalDirections(az_deg,zen_deg)*p;
assert(norm(reference)>1e-8,'FIM:AmbiguousCDLRayReference', ...
    'The CDL template direction center is ambiguous.');
reference = reference/norm(reference);
end
