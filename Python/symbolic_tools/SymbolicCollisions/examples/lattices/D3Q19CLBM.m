%--------------------------------------------------------------------------
% This program is part of the paper "Three-dimensional cascaded lattice 
% Boltzmann method: improved implementation and consistent forcing scheme".
% It allows the readers to calculate M, N, M^-1 and N^-1 for  D3Q19 CLBM
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% Copyright (C) 2018  Linlin Fei (fll15@mails.tsinghua.edu.cn)
% Address: Center for Combustion Energy; Key laboratory for Thermal Science
% and Power Engineering of Ministry of Education, Department of Energy and 
%Power Engineering, Tsinghua University, Beijing 100084
%--------------------------------------------------------------------------
clear all
clc
ex=[0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0];
ey=[0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1];
ez=[0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1];
M=zeros(19,19);
M(1,:)=ex.^0;
M(2,:)=ex;
M(3,:)=ey;
M(4,:)=ez;
M(5,:)=ex.*ey;
M(6,:)=ex.*ez;
M(7,:)=ey.*ez;
M(8,:)=ex.*ex;
M(9,:)=ey.*ey;
M(10,:)=ez.*ez;
M(11,:)=ex.*ey.*ey;
M(12,:)=ex.*ez.*ez;
M(13,:)=ey.*ex.*ex;
M(14,:)=ez.*ex.*ex;
M(15,:)=ey.*ez.*ez;
M(16,:)=ez.*ey.*ey;
M(17,:)=ex.*ex.*ey.*ey;
M(18,:)=ex.*ex.*ez.*ez;
M(19,:)=ey.*ey.*ez.*ez;
%A=M^-1;
syms ux uy uz
N(2,:)=[-ux,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
N(1,:)=[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
N(3,:)=[-uy,0,1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
N(4,:)=[-uz,0,0,1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
N(5,:)=[ux*uy,-uy,-ux,0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
N(6,:)=[ux*uz,-uz,0,-ux,0,1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0];
N(7,:)=[uy*uz,0,-uz,-uy,0,0,1,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
N(8,:)=[ux*ux,-2*ux,0,0,0,0,0,1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
N(9,:)=[uy*uy,0,-2*uy,0,0,0,0,0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
N(10,:)=[uz*uz,0,0,-2*uz,0,0,0,0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0];
N(11,:)=[-ux*uy*uy,uy*uy,2*ux*uy,0,-2*uy,0,0,0,-ux,0,1, 0, 0, 0, 0, 0, 0, 0, 0];
N(12,:)=[-ux*uz*uz,uz*uz, 0,2*ux*uz,0,-2*uz,0,0,0,-ux,0,1, 0, 0, 0, 0, 0, 0, 0];
N(13,:)=[-ux*ux*uy,2*ux*uy,ux*ux,0,-2*ux,0,0,-uy,0,0, 0,0, 1, 0, 0, 0, 0, 0, 0];
N(14,:)=[-ux*ux*uz,2*ux*uz,0, ux*ux,0,-2*ux,0,-uz,0,0, 0,0, 0, 1, 0, 0, 0, 0, 0];
N(15,:)=[-uy*uz*uz,0,uz*uz,2*uy*uz,0,0,-2*uz,0,0,-uy,0,0,0,0,1,0,0,0,0];
N(16,:)=[-uy*uy*uz,0,2*uy*uz,uy*uy,0,0,-2*uy,0,-uz,0,0,0,0,0,0,1,0,0,0];
N(17,:)=[ux*ux*uy*uy,-2*ux*uy*uy,-2*uy*ux*ux,0,4*ux*uy,0,0,uy*uy,ux*ux,0,-2*ux,0,-2*uy,0,0,0,1,0,0];
N(18,:)=[ux*ux*uz*uz,-2*ux*uz*uz,0,-2*uz*ux*ux,0,4*ux*uz,0,uz*uz,0,ux*ux,0,-2*ux,0,-2*uz,0,0,0,1,0];
N(19,:)=[uy*uy*uz*uz,0,-2*uy*uz*uz,-2*uz*uy*uy,0,0,4*uy*uz,0,uz*uz,uy*uy,0,0,0,0,-2*uy,-2*uz,0,0,1];
%B=N^-1;