function raw = generateCDLBase(modelType,varargin)
%GENERATECDLBASE Build a rawCDL template from the 3GPP CDL cluster tables.
% CenterFrequency is in Hz and DelaySpread is in ns (default 30 ns).
% Returned delays are in seconds, powers are linear, and angles are in radians.

delayProfile = upper(char(modelType));
if startsWith(delayProfile,'NR-')
    delayProfile = delayProfile(4:end);
end
delayProfile = validatestring(delayProfile,{'CDL-A','CDL-B','CDL-C','CDL-D','CDL-E'}, ...
    mfilename,'modelType');

options = inputParser;
options.FunctionName = mfilename;
options.PartialMatching = false;
addParameter(options,'CenterFrequency',3.5e9);
addParameter(options,'DelaySpread',[]);
parse(options,varargin{:});
% The tabulated normalized delays, gains and angles are frequency independent.
validateattributes(options.Results.CenterFrequency,{'numeric'}, ...
    {'real','finite','scalar','positive'},mfilename,'CenterFrequency');
delaySpread_ns = options.Results.DelaySpread;
if isempty(delaySpread_ns)
    delaySpread_ns = 30;
else
    validateattributes(delaySpread_ns,{'numeric'},{'real','finite','scalar','nonnegative'}, ...
        mfilename,'DelaySpread');
end
delaySpread_s = double(delaySpread_ns)*1e-9;

[paths,xpr_dB,firstClusterKFactor_dB] = cdlProfileTable(delayProfile);
powers = 10.^(paths(:,2).'/10);
raw.name = ['NR-' delayProfile];
raw.source = '3GPP TR 38.901 V18.0.0 CDL cluster tables';
raw.numClusters = size(paths,1);
raw.DS = delaySpread_s;
raw.KF = firstClusterKFactor_dB;
raw.XPR = xpr_dB;
raw.taus = paths(:,1).'*delaySpread_s;
raw.pow = powers/sum(powers);
raw.AoD = deg2rad(paths(:,3).');
raw.AoA = deg2rad(paths(:,4).');
raw.EoD = deg2rad(90-paths(:,5).');
raw.EoA = deg2rad(90-paths(:,6).');
end
