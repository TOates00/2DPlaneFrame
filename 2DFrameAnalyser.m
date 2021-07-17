% 2D Plane Frame 


clear 
% taking inputs from function at bottom
[nne,ndn,nns,xy,nelm,nod,emod,area,I,nnl,dnl,nbn,dbn,dml,nml] = frameinput;

% calcing number of DOF of entire structure
ntdof=nns*ndn;

% initialising global stiffness matrix, global displacement vector
gsm=zeros(ntdof,ntdof);
gdv=zeros(ntdof);

% calculating EI vector of all members (kN/mm^2 * mm^4 = kNmm^2)
EI=emod.*I;

%good. Going through each member and calcing the local stiffness matrix,,
%transformation matrix, and thus global stiffness matrix.
for ielm = 1:nelm
    % calculating lengths (in mm)
    n1=nod(ielm,1);
    n2=nod(ielm,2);
    lx=xy(n2,1)-xy(n1,1);
    ly=xy(n2,2)-xy(n1,2);
    l=sqrt(lx^2+ly^2);
    
    % Element stiffness matrix in local axis (in kN & mm)
    esml(1,1)=area(ielm)*emod(ielm)/l; 
    esml(2,1)=0;
    esml(3,1)=0;
    esml(4,1)=-esml(1,1);
    esml(5,1)=0;
    esml(6,1)=0; 
    esml(1,2)=0;
    esml(2,2)=12*EI(ielm)/l^3;
    esml(3,2)=6*EI(ielm)/l^2;
    esml(4,2)=0;
    esml(5,2)=-12*EI(ielm)/l^3;
    esml(6,2)=6*EI(ielm)/l^2; 
    esml(1,3)=0;
    esml(2,3)=6*EI(ielm)/l^2;
    esml(3,3)=4*EI(ielm)/l;
    esml(4,3)=0;
    esml(5,3)=-6*EI(ielm)/l^2;
    esml(6,3)=2*EI(ielm)/l; 
    esml(1,4)=-esml(1,1);
    esml(2,4)=0;
    esml(3,4)=0;
    esml(4,4)=esml(1,1);
    esml(5,4)=0;
    esml(6,4)=0; 
    esml(1,5)=0;
    esml(2,5)=-12*EI(ielm)/l^3;
    esml(3,5)=-6*EI(ielm)/l^2;
    esml(4,5)=0;
    esml(5,5)=12*EI(ielm)/l^3;
    esml(6,5)=-6*EI(ielm)/l^2; 
    esml(1,6)=0;
    esml(2,6)=6*EI(ielm)/l^2;
    esml(3,6)=2*EI(ielm)/l;
    esml(4,6)=0;
    esml(5,6)=-6*EI(ielm)/l^2;
    esml(6,6)=4*EI(ielm)/l; 
    
    % Transformation matrix (from local to global)
    t=zeros(6,6);          
    t(1,1)=lx/l;
    t(2,1)=-ly/l;
    t(1,2)=ly/l;
    t(2,2)=lx/l;
    t(3,3)=1;
    t(4,4)=lx/l;
    t(5,4)=-ly/l;
    t(4,5)=ly/l;
    t(5,5)=lx/l;
    t(6,6)=1;
    
    % Element stiffness matrix in global axis
    esm=t'*esml*t;     
    
    % making global stiffness matrix
    for in=1:nne
        for id=1:ndn
            il=(in-1)*ndn+id;
            ig=(nod(ielm,in)-1)*ndn+id;
    
            for jn=1:nne
                for jd=1:ndn
                  jl=(jn-1)*ndn+jd;
                  jg=(nod(ielm,jn)-1)*ndn+jd;
                  gsm(ig,jg)=gsm(ig,jg)+esm(il,jl);
                end
            end
        end
    end 
    
    
    
    
    
end

% making diagonals non-zero. Arbitrary number, so long as we get
% non-singular matrix (need non-zero diagonals)
for i=1:ntdof
    if(gsm(i,i)==0)
        gsm(i,i)=100.0;
    end
end

% checking boundary conditions, apply large stiffness coefficients at
% boundaries (restrains DOF to 0)
for ibn=1:nbn
    for id=1:ndn
        ii=(dbn(ibn,1)-1)*ndn+id;
        
        if(dbn(ibn,id+1)>0)
            gsm(ii,ii)=gsm(ii,ii)*1.0e6;
        end
    end
end

% initialising the global load matrix and local/global load vectors
glm = zeros(nns*ndn,nelm+1);
pl = zeros(ndn*nne,nelm);
pg = zeros(ndn*nne,nelm);
lm = zeros(ndn*nne,nelm);

% member load calcs
 for in1 = 1:nml
    
    % specifying the member the member load is acting in 
    ielm = dml(in1, 8);
    
    % initialising fixed end forces/moments
    femAB = 0;
    femBA = 0;
    fefyAB = 0;
    fefyBA = 0;
    fefxAB = 0;
    fefxBA = 0;
    
    n1=nod(ielm,1);
    n2=nod(ielm,2);
    lx=xy(n2,1)-xy(n1,1);
    ly=xy(n2,2)-xy(n1,2);
    l=sqrt(lx^2+ly^2);
    
    % calculating local distances load is acting away from point (useful
    % only for point load within beam forces)
    xdist = dml(in1, 6)-xy(n1,1);
    ydist = dml(in1, 7)-xy(n1,2);
    
    % creating transformation matrix (3x3) for member forces
    tf = zeros(3,3);
    tf(1,1)=lx/l;
    tf(2,1)=-ly/l;
    tf(1,2)=ly/l;
    tf(2,2)=lx/l;
    tf(3,3)=1;
    
    
    if dml(in1,5) == 1
        % conversion to local axis system if inputted as global
        lmf = tf*dml(in1,1:3)';
        coords = tf(1:2,1:2)*[xdist,ydist]';
        % calculation of a and b distances for non-symmetric loading
        b = l - coords(1);
        a = coords(1);
    elseif dml(in1,5) == 0
        % reconfiguration of local to maintain local axis system
        lmf = dml(in1,1:3)';
        coords = [xy(n1,1); xy(n1,2)]+[xdist,ydist]'; 
        b = l - coords(1);
        a = coords(1);
    end
    
    % for point load within member
    if dml(in1, 4) == 1 
        
         % kN * mm^2 * mm / mm^2 = kNmm
         femAB = -lmf(2)*b^2*a/l^2;
         femBA = lmf(2)*a^2*b/l^2;
         
         % (kNmm + kNmm - kNmm)/mm = kN
         fefyAB = (femAB+femBA-lmf(2)*b)/(l);
         fefyBA = -fefyAB-lmf(2);
         
         % mm / mm * kN = kN
         fefxAB = -b/l*lmf(1);
         fefxBA = -a/l*lmf(1);
    end
    
    % for UDL within member
    if dml(in1, 4) == 2 
         
         % kN/m/10^3 * mm^2 = kNmm
         femAB = -lmf(2)/10^3*l^2/12;
         femBA = lmf(2)/10^3*l^2/12;
         
         % kN/m/10^3 * mm = kN
         fefyAB = -lmf(2)*l/2/10^3;
         fefyBA = -lmf(2)*l/2/10^3;
         
         % mm * kN/m/10^3 = kN
         fefxAB = -l/2*lmf(1)/10^3;
         fefxBA = -l/2*lmf(1)/10^3;
    end
    
    % for point moment
    if dml(in1, 4) == 3 
        
         % kNm * 10^3 * mm * mm / mm^2 = kNmm
         femAB = -lmf(3)*10^3*b*(2*a-b)/l^2;
         femBA = -lmf(3)*10^3*a*(2*b-a)/l^2;
         
         % kNm * 10^3 * mm * mm / mm^3 = kN
         fefyAB = -6*lmf(3)*10^3*a*b/l^3;
         fefyBA = 6*lmf(3)*10^3*a*b/l^3;
         
         fefxAB = 0;
         fefxBA = 0;
    end
    
    % for non-uniform DL - triangular loading (i.e., hydrostatic loading)
    if dml(in1, 4) == 4 
         if xdist > 0
            % peak is at 'right-most' point
            % kN/m/10^3 * mm^2 = kNmm
            femAB = -lmf(2)/10^3*l^2/30;
            femBA = lmf(2)/10^3*l^2/20;
         
            % kN/m/10^3 * mm = kN
            fefyAB = -3*lmf(2)*l/20/10^3;
            fefyBA = -7*lmf(2)*l/20/10^3;
         
            % mm * kN/m/10^3 = kN
            fefxAB = -1/6*l*lmf(1)/10^3;
            fefxBA = -2/6*l*lmf(1)/10^3;
         else
            % peak is at 'left-most' point
            % kN/m/10^3 * mm^2 = kNmm
            femAB = -lmf(2)/10^3*l^2/20;
            femBA = lmf(2)/10^3*l^2/30;
         
            % kN/m/10^3 * mm = kN
            fefyAB = -7*lmf(2)*l/20/10^3;
            fefyBA = -3*lmf(2)*l/20/10^3;
         
            % mm * kN/m/10^3 = kN
            fefxAB = -2/6*l*lmf(1)/10^3;
            fefxBA = -1/6*l*lmf(1)/10^3;
         end
    end
    
    
    % creating transformation matrix (6x6) of member that forces acts in
    t=zeros(6,6);          
    t(1,1)=lx/l;
    t(2,1)=-ly/l;
    t(1,2)=ly/l;
    t(2,2)=lx/l;
    t(3,3)=1;
    t(4,4)=lx/l;
    t(5,4)=-ly/l;
    t(4,5)=ly/l;
    t(5,5)=lx/l;
    t(6,6)=1;
    
    % centralising local member forces to a single matrix (ith columns
    % represents interal forces in ith element)
    pl(:,dml(in1,8)) = [-fefxAB;
                 -fefyAB;
                 -femAB;
                 -fefxBA;
                 -fefyBA;
                 -femBA ];  
    
    % global member forces 
    pg(:,dml(in1,8)) = t'*pl(:,dml(in1,8));
    
    % updating global load matrix with global member forces for each
    % member. ii1 adds member forces for first node of member within and
    % ii2 adds member forces for second node of member within (works for 
    % structures that have more than two beams connecting at a single 
    % node).
    ii1 = 1+3*(nod(dml(in1,8),1)-1);
    ii2 = 1+3*(nod(dml(in1,8),2)-1);
    glm(ii1:(ii1+ndn-1),dml(in1,8)) = pg(1:3,dml(in1,8));
    glm(ii2:(ii2+ndn-1),dml(in1,8)) = pg(4:6,dml(in1,8));
  end

% point load calcs
for inl=1:nnl
     ii=(dnl(inl,2)-1)*ndn+dnl(inl,3);
    % adding point load calcs to global load matrix
    glm(ii,nelm+1)=glm(ii,nelm+1)+dnl(inl,1);
end

% converting global load matrix to global load vector. Adding the rows of
% the global load matrix to form the global load vector.
glv = sum(glm,2);


% Have gsm (kN & mm) * gdv(rads & mm) = glv (kN & mm)
gdv=gsm\glv;    % Nodal displacement vector of the structure

% prints nodal displacement titles
fid=fopen('frameoutput.txt','w');
fprintf(fid,'%30s\r\n','Nodal Displacement Vector');
fprintf(fid,'%16s %12s %12s %12s\n','Node Number','U','V','Theta');

% prints nodal displacement values for each member (vertical and horizontal
% displacement, and rotation)
for in=1:nns
    ii=(in-1)*ndn;
    fprintf(fid,'%16d %12.5e %12.5e %12.5e\n',in,gdv(ii+1),gdv(ii+2),gdv(ii+3));
end

% prints member end moment titles
 fprintf(fid,'%30s\r\n','Member End Forces/Moments');
 fprintf(fid,'%19s %10s %12s %12s %12s %12s %12s\n','Element Number','X12','Y12','M12','X21','Y21','M21');

% calculates local member forces and prints them 
for ielm=1:nelm
    for in=1:nne
        for id=1:ndn
            il=(in-1)*ndn+id;
            ig=(nod(ielm,in)-1)*ndn+id;
            edv(il)=gdv(ig);       % Nodal displacement vector of an element in global axis system

        end
    end
                
    n1=nod(ielm,1);
    n2=nod(ielm,2);
    lx=xy(n2,1)-xy(n1,1);
    ly=xy(n2,2)-xy(n1,2);
    l=sqrt(lx^2+ly^2);
    esml(1,1)=area(ielm)*emod(ielm)/l;  % Element stiffness matrix in local axis
    esml(2,1)=0;
    esml(3,1)=0;
    esml(4,1)=-esml(1,1);
    esml(5,1)=0;
    esml(6,1)=0; 
    esml(1,2)=0;
    esml(2,2)=12*EI(ielm)/l^3;
    esml(3,2)=6*EI(ielm)/l^2;
    esml(4,2)=0;
    esml(5,2)=-12*EI(ielm)/l^3;
    esml(6,2)=6*EI(ielm)/l^2; 
    esml(1,3)=0;
    esml(2,3)=6*EI(ielm)/l^2;
    esml(3,3)=4*EI(ielm)/l;
    esml(4,3)=0;
    esml(5,3)=-6*EI(ielm)/l^2;
    esml(6,3)=2*EI(ielm)/l; 
    esml(1,4)=-esml(1,1);
    esml(2,4)=0;
    esml(3,4)=0;
    esml(4,4)=esml(1,1);
    esml(5,4)=0;
    esml(6,4)=0; 
    esml(1,5)=0;
    esml(2,5)=-12*EI(ielm)/l^3;
    esml(3,5)=-6*EI(ielm)/l^2;
    esml(4,5)=0;
    esml(5,5)=12*EI(ielm)/l^3;
    esml(6,5)=-6*EI(ielm)/l^2; 
    esml(1,6)=0;
    esml(2,6)=6*EI(ielm)/l^2;
    esml(3,6)=2*EI(ielm)/l;
    esml(4,6)=0;
    esml(5,6)=-6*EI(ielm)/l^2;
    esml(6,6)=4*EI(ielm)/l;  
    t=zeros(6,6);          % Transformation matrix
    t(1,1)=lx/l;
    t(2,1)=-ly/l;
    t(1,2)=ly/l;
    t(2,2)=lx/l;
    t(3,3)=1;
    t(4,4)=lx/l;
    t(5,4)=-ly/l;
    t(4,5)=ly/l;
    t(5,5)=lx/l;
    t(6,6)=1;
    
    % units in mm & rads
    edvl=t*edv';       %Nodal displacement vector of an element in local axis system

    % units in kN & mm
    mef=esml*edvl;     %Member end forces
    
    % constructs local member forces (X, Y and M for each node) in all
    % members in a matrix by taking FEFs/FEMs from member end
    % forces/moments (forces due to joint displacement/rotation)
    % units in kN & mm
    lm(:,ielm) = mef - pl(:,ielm);
   
    fprintf(fid,'%17d %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n',ielm,lm(1,ielm),lm(2,ielm),lm(3,ielm),lm(4,ielm),lm(5,ielm),lm(6,ielm));
    
end


function [nne,ndn,nns,xy,nelm,nod,emod,area,I,nnl,dnl,nbn,dbn,dml,nml] = frameinput
nne=2;  % Number of nodes per element 
ndn=3;  % Number of DOF per node
nns=5;  % Number of nodes of the structure
xy=[0 0;  % xy is the coordinate of all nodes of the structure  
    0 8000; 
    6000 12500;
    12000 8000;
    12000 0];
I=[395.7*10^6 395.7*10^6 395.7*10^6 395.7*10^6];
nelm=4;  % Number of elements
nod=[1 2; % Nodal connection matrix 
    2 3;
    3 4;
    4 5];
emod = [200 200 200 200];  % Elastic modulus of the elements
area = [15600 15600 15600 15600];  % Cross-sectional area of the elements
nnl=2;  % Number of point loads acting at the nodes directly 
nml=2;  % Number of loads acting within members

dnl =[40.0 2 1;
      -50 3 2];  % Details of the point loads (magnitude, node, direction (1 = horizontal, 2 = vertical, 3 = moment)

% Details of in-member forces: P_x, P_y, P_Mz, load type 
% (1 = point, 2 = UDL, 3 = point moment, 4 = triangular DL), axis systems 
% (local = 0, global = 1), x1, y1, member within  
dml = [2.4, -3.2, 0, 2, 1, 3000, 4000, 2;
       -50, 0, 0, 1, 1, 8400, 10700, 3];
   
% Note: axis system is for both forces and coordinates of force    
% Note: (x1,y1) is the location of the forces (or the peak force in the
% case of a non-uniform DL)

nbn=2;  % Number of boundary nodes 
dbn=[1 1 1 1;   % Detail of the boundary nodes (node, DOFs restrained)
     5 1 1 1];
end

