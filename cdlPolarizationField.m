function field = cdlPolarizationField(directions,rotation,slant_deg)
%CDLPOLARIZATIONFIELD Unit element field in reference-frame theta/phi bases.
% rotation maps element-local vectors to the frame used by the 3-by-L
% unit directions. The element power pattern is applied by the caller.
% TR 38.901 V18.0.0, clauses 7.1.2 and 7.3.2 (Model-2).

[thetaLocal,phiLocal] = sphericalBasis(rotation.'*directions);
cartesian = rotation*(cosd(slant_deg)*thetaLocal+sind(slant_deg)*phiLocal);
[theta,phi] = sphericalBasis(directions);
field = [sum(theta.*cartesian,1);sum(phi.*cartesian,1)];
end

function [theta,phi] = sphericalBasis(direction)
azimuth = atan2(direction(2,:),direction(1,:));
x = cos(azimuth);
y = sin(azimuth);
z = direction(3,:);
theta = [z.*x;z.*y;-hypot(direction(1,:),direction(2,:))];
phi = [-y;x;zeros(size(z))];
end
