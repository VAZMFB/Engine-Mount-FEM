% FEM engine mount calculation
% Author: Milos D. Petrasinovic <mpetrasinovic@mas.bg.ac.rs>
% Structural Analysis of Flying Vehicles
% Faculty of Mechanical Engineering, University of Belgrade
% Department of Aerospace Engineering, Flying structures
% https://vazmfb.com
% Belgrade, 2021
%
% ---------------
%
% Copyright (C) 2021 Milos Petrasinovic <info@vazmfb.com>
%  
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as 
% published by the Free Software Foundation, either version 3 of the 
% License, or (at your option) any later version.
%   
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%   
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% ---------------
clear('all'), clc, close('all'), tic
disp([' --- ' mfilename ' --- ']);

% - Input parameters

% Characteristics of cross section
Ds = 16; % [mm] outer diameter of the pipe
t = 1; % [mm] Pipe wall thickness
% Material characteristics
E = 210e9; % [Pa] Elastic modulus

% Coordinates of nodes
% nc(node, 1:3) = [x, y, z];
nc = [0, 0, 0; % node 1
      0, 700, 0; % node 2
      0, 0, 700; % node 3
      0, 700, 700; % node 4
      1200, 0, 350; % node 5
      1200, 700, 350; % node 6
      300, 37.5, 350; % node 7
      300, 662.5, 350; % node 8
      900, 350, 37.5; % node 9
      900, 350, 662.5;]; % node 10

% Nodes of engine mount rods
% enn(element, :) = [node1, node2];
enn = [1, 7; % rod 1
      3, 7; % rod 2
      2, 8; % rod 3
      4, 8; % rod 4
      5, 10; % rod 5
      5, 9; % rod 6
      6, 10; % rod 7
      6, 9;];  % rod 8

% Nodes of rods that simulates engine (infinitely rigid rods)
% enm(element, :) = [node1, node2];
enm = [7, 8; % engine rod 1
       8, 9; % engine rod 2
       9, 7; % engine rod 3
       7, 10; % engine rod 4
       10, 8; % engine rod 5
       10, 9;]; % engine rod 6

% Boundary conditions
% bc(node, 1:3) = [Tx, Ty, Tz];
% 0 - free
% 1 - constrained
bc(1, :) = [1, 1, 1]; % node 1
bc(2, :) = [1, 1, 1]; % node 2
bc(3, :) = [1, 1, 1]; % node 3
bc(4, :) = [1, 1, 1]; % node 4
bc(5, :) = [1, 1, 1]; % node 5
bc(6, :) = [1, 1, 1]; % node 6

% External loads
% F(node, 1:3) = [Fx[N], Fy[N], Fz[N]];
F(7, :) = [3906.95, 0, -5036.15]; % node 7
F(8, :) = [3906.95, 0, -5036.15]; % node 8
F(9, :) = [3906.95, 0, -5036.15]; % node 9
F(10, :) = [3906.95, 0, -5036.15]; % node 10

% Additional variables
s = [1000, % [m] to [mm]
     200, % displacement magnification
     100, % vector length for boundary conditions
     1/50, % magnification of the vector for external loads
     100; % free space around the model
     0.2, % the size of the vector arrow
     -60]; % moving node numbers

% - Defining a stiffness matrix
d2r = 180/pi; % degrees to radians
en = [enn; enm]; % nodes of all elements
ne = size(en, 1); % number of elements
nen = size(enn, 1);
nem = size(enm, 1);
nn = size(nc, 1); % number of nodes
dof = 3*nn; % number of degrees of freedom
D = zeros(dof, 1); % displacements
bc(end+1:nn, :) = zeros(length(size(bc, 1)+1:nn), 3); % boundary conditions
F(end+1:nn, :) = zeros(length(size(F, 1)+1:nn), 3); % external loads
F = reshape(F.', [], 1); 
K = zeros(dof, dof); % stiffness matrix
sigma = zeros(nen, 1);

% Degrees of freedom of nodes
rdof = find(bc')'; % constrained degrees of freedom
fdof = find(~bc')'; % free degrees of freedom

% Cross-sectional characteristics (circular tube)
Ds = Ds/s(1); % [m]
t = t/s(1); % [m]
A = (Ds^2-(Ds-2*t)^2)*pi/4*ones(ne, 1); % [m^2] Cross-sectional area
A(nen+1:end) = A(nen+1:end)*10^6; % increasing stiffness

% Determination of stiffness matrix
nc = nc/s(1); % [m]
for i=1:ne 
  j = en(i, :);       
  edof = [3*j(1)-2, 3*j(1)-1, 3*j(1) ... % degrees of freedom of the element
         3*j(2)-2, 3*j(2)-1, 3*j(2)]; 
  L_e = sqrt((nc(j(2), 1)-nc(j(1), 1))*(nc(j(2), 1)-... % element length
      nc(j(1), 1))+(nc(j(2), 2)-nc(j(1), 2))*(nc(j(2), 2)-...
      nc(j(1), 2))+(nc(j(2), 3)-nc(j(1), 3))*(nc(j(2), 3)-nc(j(1), 3)));
  CXx = (nc(j(2), 1)-nc(j(1), 1))/L_e;
  CYx = (nc(j(2), 2)-nc(j(1), 2))/L_e;
  CZx = (nc(j(2), 3)-nc(j(1), 3))/L_e;
  r = [CXx*CXx, CXx*CYx, CXx*CZx;
       CYx*CXx, CYx*CYx, CYx*CZx;
       CZx*CXx, CZx*CYx, CZx*CZx];
  
  % Element stiffness matrix in the local coordinate system of the element
  k_e = (E*A(i))/L_e;
  
  T = [r, -r; -r, r];  % Transformation matrix   
  
  K(edof, edof) = K(edof, edof)+k_e*T; % Element stiffness matrix
end  

% - Solving FEM
% Solving the reduced finite element equation
D1 = K(fdof, fdof)\F(fdof);
D(fdof) = D1;
F1 = K*D;
R = F1(rdof); % bond reaction vector

for i=1:nen
  j = en(i, :);       
  edof = [3*j(1)-2, 3*j(1)-1, 3*j(1) ... % degrees of freedom of the element
         3*j(2)-2, 3*j(2)-1, 3*j(2)]; 
  L_e = sqrt((nc(j(2), 1)-nc(j(1), 1))*(nc(j(2), 1)-... % element length
      nc(j(1), 1))+(nc(j(2), 2)-nc(j(1), 2))*(nc(j(2), 2)-...
      nc(j(1), 2))+(nc(j(2), 3)-nc(j(1), 3))*(nc(j(2), 3)-nc(j(1), 3)));
  CXx = (nc(j(2), 1)-nc(j(1), 1))/L_e;
  CYx = (nc(j(2), 2)-nc(j(1), 2))/L_e;
  CZx = (nc(j(2), 3)-nc(j(1), 3))/L_e;
  sigma(i) = E/L_e*[-CXx, -CYx, -CZx, CXx, CYx, CZx]*D(edof); % [Pa] Normal stress
end 

% Determination of forces in rods
Q = sigma.*A(1:nen); % [N] Foreces in rods

% - Display of results
% Display of displacement values and reactions
disp(' -------------------- ');
disp(' Displacements')
i = 1:dof;
Dn = reshape(D, 3, []); % nodes displacements
Dv = [reshape(repmat(1:nn, 3, 1), 1, []);...
  reshape(repmat(1:3, nn, 1).', [], 1).'; Dn(i)*s(1)];
disp(' Node | Component | Displacement [mm]');
fprintf(' %3d | %3d | %14.10f\n', Dv);

% Display of reactions
disp(' -------------------- ');
disp(' Reactions')
nrdof = mod(rdof, 3); % constrained degrees of freedom
nrdof(nrdof == 0) = ones(length(find(nrdof == 0)), 1)*3;
nr = (rdof-nrdof)/3+1; % node with constraints
Fv = [nr.', nrdof.', R].';
disp(' Node | Component | Reaction [N]');
fprintf(' %3d | %3d | %14.10f\n', Fv);

% Display of forces in rods
disp(' -------------------- ');
disp(' Forces in rods')
Qv = [enn(:, 1), enn(:, 2), Q].';
disp(' Rod | Force [N]');
fprintf(' %d-%d | %14.10f\n', Qv);

% - Display of model
% Display of initial model with loads
disp(' -------------------- ');
disp(' Display of initial model with loads... ');
drawArrow = @(x, y, z, varargin) quiver3(x(1), y(1), z(1), ...
    x(2)+10^-5, y(2)+10^-5, z(2)+10^-5, 0, varargin{:});   
drawArrowMarker = @(x, y, z, m, varargin) [plot3([x(1); x(1)+x(2)], ...
  [y(1); y(1)+y(2)], [z(1); z(1)+z(2)], '-', varargin{:}), ...
  plot3(x(1)+x(2), y(1)+y(2), z(1)+z(2), m, varargin{:})];    

nc = nc*s(1);
rds = [[nc(en(:, 1), 1), nc(en(:, 2), 1)], ...
    [nc(en(:, 1), 2), nc(en(:, 2), 2)], ...
    [nc(en(:, 1), 3), nc(en(:, 2), 3)]]; % rods
  
figure(1);
box on, grid on, hold on
boje = get(gca, 'colororder');
plot3(rds(1:nen, 1:2).', rds(1:nen, 3:4).', rds(1:nen, 5:6).', ...
    'LineWidth', 2, 'Color', boje(2, :));
plot3(rds(nen+1:end, 1:2).', rds(nen+1:end, 3:4).', rds(nen+1:end, 5:6).', ...
    'LineWidth', 2, 'Color', boje(6, :));
plot3(nc(:, 1), nc(:, 2), nc(:, 3), 'o', 'LineWidth', 2, ...
    'MarkerFaceColor', boje(2, :), 'Color', boje(2, :));
rotate3d, view(45, 15), axis equal, set(gca, 'FontSize', 18);
Xlims = [min(nc(:, 1)), max(nc(:, 1))];
Ylims = [min(nc(:, 2)), max(nc(:, 2))];
Zlims = [min(nc(:, 3)), max(nc(:, 3))];
xlim([Xlims(1)-s(5), Xlims(2)+s(5)]);
ylim([Ylims(1)-s(5), Ylims(2)+s(5)]);
zlim([Zlims(1)-s(5), Zlims(2)+s(5)]);

% Display of boundary conditions
for i=1:length(rdof)
  p = zeros(1, 3);
  p(nrdof(i)) = 1*s(3);
  drawArrowMarker([nc(nr(i), 1), p(1)], [nc(nr(i), 2), p(2)], ...
    [nc(nr(i), 3), p(3)], '+', 'Color', 'b', 'LineWidth', 2, ...
    'MarkerFaceColor', 'b');
end

% Display of external forces
for i=1:nn
  p = zeros(1, 3);
  if(norm(F((i-1)*3+1:(i-1)*3+3)) > 0)
    p = (F((i-1)*3+1:(i-1)*3+3))*s(4);
    drawArrow([nc(i, 1), p(1)], [nc(i, 2), p(2)], ...
      [nc(i, 3), p(3)], 'Color', 'b', 'LineWidth', 2, ...
      'MarkerFaceColor', 'b', 'MaxHeadSize', s(6));
  end
end

% Nodes numbering
text(nc(:, 1)+s(7), nc(:, 2)+s(7), nc(:, 3), ...
   num2str((1:nn).'), 'FontSize', 18, 'Color', boje(2, :));
print('fig-1.png', '-dpng', '-F:18');
print('fig-1.svg', '-dsvg', '-FCMU Serif:18');

% Display of deformed model with loads
disp(' Display of deformed model with loads... ');
ncd = nc+(Dn.'.*s(1)).*s(2);
rdsd = [[ncd(en(:, 1), 1), ncd(en(:, 2), 1)], ...
    [ncd(en(:, 1), 2), ncd(en(:, 2), 2)], ...
    [ncd(en(:, 1), 3), ncd(en(:, 2), 3)]]; % rods

figure(2);
box on, grid on, hold on
boje = get(gca, 'colororder');
plot3(rds(:, 1:2).', rds(:, 3:4).', rds(:, 5:6).', '--', ...
    'LineWidth', 1, 'Color', boje(1, :));
plot3(rdsd(1:nen, 1:2).', rdsd(1:nen, 3:4).', rdsd(1:nen, 5:6).', ...
    'LineWidth', 2, 'Color', boje(2, :));
plot3(rdsd(nen+1:end, 1:2).', rdsd(nen+1:end, 3:4).', rdsd(nen+1:end, 5:6).', ...
    'LineWidth', 2, 'Color', boje(6, :));
plot3(ncd(:, 1), ncd(:, 2), ncd(:, 3), 'o', 'LineWidth', 2, ...
    'MarkerFaceColor', boje(2, :), 'Color', boje(2, :));
rotate3d, view(45, 15), axis equal, set(gca, 'FontSize', 18);
xlim([Xlims(1)-s(5), Xlims(2)+s(5)]);
ylim([Ylims(1)-s(5), Ylims(2)+s(5)]);
zlim([Zlims(1)-s(5), Zlims(2)+s(5)]);
print('fig-2.png', '-dpng', '-F:18');
print('fig-2.svg', '-dsvg', '-FCMU Serif:18');

% - End of program
disp(' -------------------- ');
disp(' The program was successfully executed... ');
disp([' Execution time: ' num2str(toc, '%.2f') ' seconds']);
disp(' -------------------- ');