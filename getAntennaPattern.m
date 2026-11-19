function [Vgain,Vphase] = getAntennaPattern(Phi,Theta)
%GETANTENNAPATTERN 保留用户接口的笛卡尔积行为，输出 numel(Theta) x numel(Phi)。
% 按原实现：Phi 0:359 先变化，Theta 0:180 后变化。不是一一配对查询。
persistent Phi_grid Theta_grid Vgain_grid Vphase_grid
if isempty(Vgain_grid)
    dataFile = fullfile(fileparts(mfilename('fullpath')),'HW_Airview_Fre6700M_E06_1to6_12.7dB.txt');
    fid = fopen(dataFile,'r');
    if fid==-1
        error('Cannot open measured pattern: %s',dataFile);
    end
    cleanup = onCleanup(@() fclose(fid));
    fgetl(fid);
    fgetl(fid);
    data = fscanf(fid,'%f %f %f %f',[4 inf]).';
    assert(size(data,1)==360*181,'Expected 360*181 pattern rows after two header lines.');
    [P,T] = ndgrid(0:359,0:180);
    assert(all(abs(data(:,1)-P(:))<1e-8) && all(abs(data(:,2)-T(:))<1e-8),'Pattern order must be Phi-fast, Theta-slow.');
    Vgain_grid = reshape(data(:,3),360,181);
    Vphase_grid = reshape(data(:,4),360,181);
    Phi_grid = 0:359;
    Theta_grid = 0:180;
end
Phi = reshape(Phi,1,[]);
Theta = reshape(Theta,1,[]);
[Phi_mesh,Theta_mesh] = meshgrid(Phi,Theta);
Vgain = interp2(Theta_grid,Phi_grid,Vgain_grid,Theta_mesh,Phi_mesh,'linear',NaN);
Vphase = interp2(Theta_grid,Phi_grid,Vphase_grid,Theta_mesh,Phi_mesh,'linear',NaN);
if any(Phi<0 | Phi>=360) || any(Theta<0 | Theta>180)
    warning('Some queries are outside Phi=[0,360), Theta=[0,180].');
end
% 原接口的周期接缝和 dB 插值行为保留；后续旋转查询由独立模块在功率域处理。
end
