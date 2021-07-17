function [nne,ndn,nns,xy,nelm,nod,emod,area,I,nnl,dnl,nbn,dbn,dml,nml] = frameinput2021
nne=2;  % Number of nodes per element 
ndn=3;  % Number of DOF per node
nns=4;  % Number of nodes of the structure 
xy=[0 0;  % xy is the coordinate of all nodes of the structure  
    4500 6000; 
    8500 6000;
    16500 0];
I=[2.923*10^7 2.923*10^7 2.923*10^7];
nelm=3;  % Number of elements
nod=[1 2; % Nodal connection matrix 
    2 3;
    3 4];
emod = [200 200 200];  % Elastic modulus of the elements
area = [4656 4656 4656];  % Cross-sectional area of the elements
nnl=1;  % Number of point loads acting at the nodes directly 
nml=3;  % Number of loads acting within members

% Details of the point loads: Magnitude, Node, Direction (1 = horizontal, 
% 2 = vertical, 3 = moment)
dnl =[20.0 2 1];  

% Details of in-member forces: P_x1, P_y1, P_Mz1, P_x2, P_y2, P_Mz2, load type 
% (1 = point, 2 = UDL, 3 = point moment, 4 = triangular DL), axis systems 
% (local = 0, global = 1), x1, y1, member within
dml = [15, 0, 0, 1, 1, 3000, 4000, 1;
        0, -4, 0, 2, 1, 4500, 6000, 2;
        0, -2.25, 0, 2, 1, 8500, 6000, 3];

% Note: axis system is for both forces and coordinates of force    
% Note: (x1,y1) is the location of the forces (or the peak force in the
% case of a non-uniform DL)

% Loads (P_X, P_Y, P_MZ), load type (1 = point, 2 = UDL, 3 = point moment), direction
% (degrees), x1, y1, x2, y2, member within


    % Details for member loads
nbn=2;  % Number of boundary nodes 
dbn=[1 1 1 1;   % Detail of the boundary nodes (node, DOFs restrained)
     4 1 1 1];
end