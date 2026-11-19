function [field,initialPhases,doppler_Hz] = cdlRayCoefficients(rays,txDirections,txRotation, ...
    txSlant_deg,velocityTx_mps,time_s,carrier_Hz,lightSpeed_mps,losDistance_m,initialPhases)
%CDLRAYCOEFFICIENTS One transmitted port's incident fields before reception.
% The 2-by-L result is the zenith/azimuth field in each ray's arrival basis
% in the same reference frame as txDirections, txRotation and velocityTx.
% initialPhases rows are theta-theta, theta-phi, phi-theta and phi-phi.
% Pass initialPhases back to evaluate a different time on the same link.
% doppler_Hz is a row vector with one Doppler frequency per ray.
% TR 38.901 V18.0.0, clauses 7.5 Steps 10-11 and 7.6.10 (mobile Tx).

los = rays.isLOS(:).';
if nargin<10 || isempty(initialPhases)
    initialPhases = zeros(4,numel(los));
    initialPhases(:,~los) = 2*pi*rand(4,nnz(~los))-pi;
else
    validateattributes(initialPhases,{'numeric'}, ...
        {'real','finite','size',[4,numel(los)]},mfilename,'initialPhases');
end

transmitted = cdlPolarizationField(txDirections,txRotation,txSlant_deg);
phase = exp(1i*initialPhases);
crossAmplitude = 1./sqrt(rays.xpr(:).');
field = [
    phase(1,:).*transmitted(1,:) + crossAmplitude.*phase(2,:).*transmitted(2,:);
    crossAmplitude.*phase(3,:).*transmitted(1,:) + phase(4,:).*transmitted(2,:);
];

carrierPerMeter = carrier_Hz/lightSpeed_mps;
if any(los)
    losPhase = exp(-1i*2*pi*carrierPerMeter*losDistance_m);
    field(:,los) = [1;-1].*transmitted(:,los)*losPhase;
end
doppler_Hz = carrierPerMeter*(velocityTx_mps(:).'*txDirections);
field = field.*(sqrt(rays.power(:).').*exp(1i*2*pi*time_s*doppler_Hz));
end
