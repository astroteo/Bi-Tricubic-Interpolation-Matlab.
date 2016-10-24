function [Phi_int,f,H]= Interpl_3D(X,X_max,Y_max,Z_max,step_grid,Phi_3d,Interp_type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%X: 3D point at wich the interpolation is computed
%{X_max,Y_max,Z_max}: extremes coordinates of the grid (Half lenght of the box).

%step_grid: dimension of the elements of  the equipaced grid.
%  CHOOSE: Interp_type=0 ---> nearest neighbor interpolation
%          Interp_type=1 ---> tri-linear interpolation (function will
%                                                    belong to C-1 class.)
%          Interp_type=2 -->tri-cubic interpolation (function will belong
%                                                           to C-2 class.)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x=X(1);
y=X(2);
z=X(3);

R1_boundedx=[-X_max:step_grid:X_max];
R1_boundedy=[-Y_max:step_grid:Y_max];
R1_boundedz=[-Z_max:step_grid:Z_max];

i_bottom=floor((x+X_max)/step_grid);
i_top=i_bottom+1;
        
        
j_bottom=floor((y+Y_max)/step_grid);
j_top=j_bottom+1;

k_bottom=floor((z+Z_max)/step_grid);
k_top=k_bottom+1;
      
x_bottom=R1_boundedx(i_bottom);
y_bottom=R1_boundedy(j_bottom);
z_bottom=R1_boundedz(k_bottom);



if Interp_type==2 % bicubic interpolation

% function values at corners: [8 Values]   

f_000=Phi_3d(i_bottom,j_bottom,k_bottom);
f_100=Phi_3d(i_top,j_bottom,k_bottom);
f_010=Phi_3d(i_bottom,j_top,k_bottom);
f_001=Phi_3d(i_bottom,j_bottom,k_top);
f_110=Phi_3d(i_top,j_top,k_bottom);
f_101=Phi_3d(i_top,j_bottom,k_top);
f_011=Phi_3d(i_bottom,j_top,k_top);
f_111=Phi_3d(i_top,j_top,k_top);

% derivative values at corners: [8x3=24 Values]

fx_000=1/(2*step_grid)*(Phi_3d(i_top,j_bottom,k_bottom)-Phi_3d(i_bottom-1,j_bottom,k_bottom));
fy_000=1/(2*step_grid)*(Phi_3d(i_bottom,j_top,k_bottom)-Phi_3d(i_bottom,j_bottom-1,k_bottom));
fz_000=1/(2*step_grid)*(Phi_3d(i_bottom,j_bottom,k_top)-Phi_3d(i_bottom,j_bottom,k_bottom-1));

fx_100=1/(2*step_grid)*(Phi_3d(i_top+1,j_bottom,k_bottom)-Phi_3d(i_bottom,j_bottom,k_bottom));
fy_100=1/(2*step_grid)*(Phi_3d(i_top,j_top,k_bottom)-Phi_3d(i_top,j_bottom-1,k_bottom));
fz_100=1/(2*step_grid)*(Phi_3d(i_top,j_bottom,k_top)-Phi_3d(i_top,j_bottom,k_bottom-1));

fx_010=1/(2*step_grid)*(Phi_3d(i_top,j_top,k_bottom)-Phi_3d(i_bottom-1,j_top,k_bottom));
fy_010=1/(2*step_grid)*(Phi_3d(i_bottom,j_top+1,k_bottom)-Phi_3d(i_bottom,j_bottom,k_bottom));
fz_010=1/(2*step_grid)*(Phi_3d(i_bottom,j_top,k_top)-Phi_3d(i_bottom,j_top,k_bottom-1));

fx_001=1/(2*step_grid)*(Phi_3d(i_top,j_bottom,k_top)-Phi_3d(i_bottom-1,j_bottom,k_top));
fy_001=1/(2*step_grid)*(Phi_3d(i_bottom,j_top,k_top)-Phi_3d(i_bottom,j_bottom-1,k_top));
fz_001=1/(2*step_grid)*(Phi_3d(i_bottom,j_bottom,k_top+1)-Phi_3d(i_bottom,j_bottom,k_bottom));

fx_110=1/(2*step_grid)*(Phi_3d(i_top+1,j_top,k_bottom)-Phi_3d(i_bottom,j_top,k_bottom));
fy_110=1/(2*step_grid)*(Phi_3d(i_top,j_top+1,k_bottom)-Phi_3d(i_top,j_bottom,k_bottom));
fz_110=1/(2*step_grid)*(Phi_3d(i_top,j_top,k_top)-Phi_3d(i_top,j_top,k_bottom-1));

fx_101=1/(2*step_grid)*(Phi_3d(i_top+1,j_bottom,k_top)-Phi_3d(i_bottom,j_bottom,k_top));
fy_101=1/(2*step_grid)*(Phi_3d(i_top,j_top,k_top)-Phi_3d(i_top,j_bottom-1,k_top));
fz_101=1/(2*step_grid)*(Phi_3d(i_top,j_bottom,k_top+1)-Phi_3d(i_top,j_bottom,k_bottom));

fx_011=1/(2*step_grid)*(Phi_3d(i_top,j_top,k_top)-Phi_3d(i_bottom-1,j_top,k_top));
fy_011=1/(2*step_grid)*(Phi_3d(i_bottom,j_top+1,k_top)-Phi_3d(i_bottom,j_bottom,k_top));
fz_011=1/(2*step_grid)*(Phi_3d(i_bottom,j_top,k_top+1)-Phi_3d(i_bottom,j_top,k_bottom));

fx_111=1/(2*step_grid)*(Phi_3d(i_top+1,j_top,k_top)-Phi_3d(i_bottom,j_top,k_top));
fy_111=1/(2*step_grid)*(Phi_3d(i_top,j_top+1,k_top)-Phi_3d(i_top,j_bottom,k_top));
fz_111=1/(2*step_grid)*(Phi_3d(i_top,j_top,k_top+1)-Phi_3d(i_top,j_top,k_bottom));
 
% cross derivatives at corners:[3x8=24 + 8x1=8 == 32 values]

%{0 0 0}
i_bb=i_bottom;
j_bb=j_bottom;
k_bb=k_bottom;

i_tt=i_top;
j_tt=j_top;
k_tt=k_top;

fxy_000=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb-1,k_bb)- Phi_3d(i_bb,j_bb-1,k_bb) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)+ Phi_3d(i_tt,j_tt,k_bb));
fxz_000=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_tt,j_bb,k_tt));
fyz_000=1/(2*step_grid^2)*( Phi_3d(i_bb,j_bb-1,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb,j_bb-1,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_bb,j_tt,k_tt));

fxyz_000=1/(8*step_grid^3)*(Phi_3d(i_tt,j_tt,k_tt)-Phi_3d(i_tt,j_bb-1,k_tt)-Phi_3d(i_bb-1,j_tt,k_tt)+Phi_3d(i_bb-1,j_bb-1,k_tt)-Phi_3d(i_tt,j_tt,k_bb-1)+Phi_3d(i_tt,j_bb-1,k_bb-1)+Phi_3d(i_bb-1,j_tt,k_bb-1)-Phi_3d(i_bb-1,j_bb-1,k_bb-1));

%{1 0 0}
i_bb=i_top;
j_bb=j_bottom;
k_bb=k_bottom;

i_tt=i_top+1;
j_tt=j_top;
k_tt=k_top;

fxy_100=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb-1,k_bb)- Phi_3d(i_bb,j_bb-1,k_bb) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)+ Phi_3d(i_tt,j_tt,k_bb));
fxz_100=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_tt,j_bb,k_tt));
fyz_100=1/(2*step_grid^2)*( Phi_3d(i_bb,j_bb-1,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb,j_bb-1,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_bb,j_tt,k_tt));

fxyz_100=1/(8*step_grid^3)*(Phi_3d(i_tt,j_tt,k_tt)-Phi_3d(i_tt,j_bb-1,k_tt)-Phi_3d(i_bb-1,j_tt,k_tt)+Phi_3d(i_bb-1,j_bb-1,k_tt)-Phi_3d(i_tt,j_tt,k_bb-1)+Phi_3d(i_tt,j_bb-1,k_bb-1)+Phi_3d(i_bb-1,j_tt,k_bb-1)-Phi_3d(i_bb-1,j_bb-1,k_bb-1));

%{0 1 0}
i_bb=i_bottom;
j_bb=j_top;
k_bb=k_bottom;

i_tt=i_top;
j_tt=j_top+1;
k_tt=k_top;

fxy_010=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb-1,k_bb)- Phi_3d(i_bb,j_bb-1,k_bb) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)+ Phi_3d(i_tt,j_tt,k_bb));
fxz_010=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_tt,j_bb,k_tt));
fyz_010=1/(2*step_grid^2)*( Phi_3d(i_bb,j_bb-1,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb,j_bb-1,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_bb,j_tt,k_tt));

fxyz_010=1/(8*step_grid^3)*(Phi_3d(i_tt,j_tt,k_tt)-Phi_3d(i_tt,j_bb-1,k_tt)-Phi_3d(i_bb-1,j_tt,k_tt)+Phi_3d(i_bb-1,j_bb-1,k_tt)-Phi_3d(i_tt,j_tt,k_bb-1)+Phi_3d(i_tt,j_bb-1,k_bb-1)+Phi_3d(i_bb-1,j_tt,k_bb-1)-Phi_3d(i_bb-1,j_bb-1,k_bb-1));

%{0 0 1}
i_bb=i_bottom;
j_bb=j_bottom;
k_bb=k_top;

i_tt=i_top;
j_tt=j_top;
k_tt=k_top+1;

fxy_001=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb-1,k_bb)- Phi_3d(i_bb,j_bb-1,k_bb) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)+ Phi_3d(i_tt,j_tt,k_bb));
fxz_001=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_tt,j_bb,k_tt));
fyz_001=1/(2*step_grid^2)*( Phi_3d(i_bb,j_bb-1,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb,j_bb-1,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_bb,j_tt,k_tt));


fxyz_001=1/(8*step_grid^3)*(Phi_3d(i_tt,j_tt,k_tt)-Phi_3d(i_tt,j_bb-1,k_tt)-Phi_3d(i_bb-1,j_tt,k_tt)+Phi_3d(i_bb-1,j_bb-1,k_tt)-Phi_3d(i_tt,j_tt,k_bb-1)+Phi_3d(i_tt,j_bb-1,k_bb-1)+Phi_3d(i_bb-1,j_tt,k_bb-1)-Phi_3d(i_bb-1,j_bb-1,k_bb-1));

%{1 1 0}
i_bb=i_top;
j_bb=j_top;
k_bb=k_bottom;

i_tt=i_top+1;
j_tt=j_top+1;
k_tt=k_top;

fxy_110=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb-1,k_bb)- Phi_3d(i_bb,j_bb-1,k_bb) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)+ Phi_3d(i_tt,j_tt,k_bb));
fxz_110=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_tt,j_bb,k_tt));
fyz_110=1/(2*step_grid^2)*( Phi_3d(i_bb,j_bb-1,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb,j_bb-1,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_bb,j_tt,k_tt));

fxyz_110=1/(8*step_grid^3)*(Phi_3d(i_tt,j_tt,k_tt)-Phi_3d(i_tt,j_bb-1,k_tt)-Phi_3d(i_bb-1,j_tt,k_tt)+Phi_3d(i_bb-1,j_bb-1,k_tt)-Phi_3d(i_tt,j_tt,k_bb-1)+Phi_3d(i_tt,j_bb-1,k_bb-1)+Phi_3d(i_bb-1,j_tt,k_bb-1)-Phi_3d(i_bb-1,j_bb-1,k_bb-1));

%{1 0 1}
i_bb=i_top;
j_bb=j_bottom;
k_bb=k_top;

i_tt=i_top+1;
j_tt=j_top;
k_tt=k_top+1;

fxy_101=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb-1,k_bb)- Phi_3d(i_bb,j_bb-1,k_bb) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)+ Phi_3d(i_tt,j_tt,k_bb));
fxz_101=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_tt,j_bb,k_tt));
fyz_101=1/(2*step_grid^2)*( Phi_3d(i_bb,j_bb-1,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb,j_bb-1,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_bb,j_tt,k_tt));

fxyz_101=1/(8*step_grid^3)*(Phi_3d(i_tt,j_tt,k_tt)-Phi_3d(i_tt,j_bb-1,k_tt)-Phi_3d(i_bb-1,j_tt,k_tt)+Phi_3d(i_bb-1,j_bb-1,k_tt)-Phi_3d(i_tt,j_tt,k_bb-1)+Phi_3d(i_tt,j_bb-1,k_bb-1)+Phi_3d(i_bb-1,j_tt,k_bb-1)-Phi_3d(i_bb-1,j_bb-1,k_bb-1));

%{0 1 1}
i_bb=i_bottom;
j_bb=j_top;
k_bb=k_top;

i_tt=i_top;
j_tt=j_top+1;
k_tt=k_top+1;

fxy_011=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb-1,k_bb)- Phi_3d(i_bb,j_bb-1,k_bb) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)+ Phi_3d(i_tt,j_tt,k_bb));
fxz_011=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_tt,j_bb,k_tt));
fyz_011=1/(2*step_grid^2)*( Phi_3d(i_bb,j_bb-1,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb,j_bb-1,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_bb,j_tt,k_tt));


fxyz_011=1/(8*step_grid^3)*(Phi_3d(i_tt,j_tt,k_tt)-Phi_3d(i_tt,j_bb-1,k_tt)-Phi_3d(i_bb-1,j_tt,k_tt)+Phi_3d(i_bb-1,j_bb-1,k_tt)-Phi_3d(i_tt,j_tt,k_bb-1)+Phi_3d(i_tt,j_bb-1,k_bb-1)+Phi_3d(i_bb-1,j_tt,k_bb-1)-Phi_3d(i_bb-1,j_bb-1,k_bb-1));

%{1 1 1}
i_bb=i_top;
j_bb=j_top;
k_bb=k_top;

i_tt=i_top+1;
j_tt=j_top+1;
k_tt=k_top+1;

fxy_111=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb-1,k_bb)- Phi_3d(i_bb,j_bb-1,k_bb) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)+ Phi_3d(i_tt,j_tt,k_bb));
fxz_111=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb-1,j_bb,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_tt,j_bb,k_tt));
fyz_111=1/(2*step_grid^2)*( Phi_3d(i_bb,j_bb-1,k_bb-1)- Phi_3d(i_bb,j_bb,k_bb-1) - Phi_3d(i_bb,j_bb-1,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_bb,j_tt,k_tt));

fxyz_111=1/(8*step_grid^3)*(Phi_3d(i_tt,j_tt,k_tt)-Phi_3d(i_tt,j_bb-1,k_tt)-Phi_3d(i_bb-1,j_tt,k_tt)+Phi_3d(i_bb-1,j_bb-1,k_tt)-Phi_3d(i_tt,j_tt,k_bb-1)+Phi_3d(i_tt,j_bb-1,k_bb-1)+Phi_3d(i_bb-1,j_tt,k_bb-1)-Phi_3d(i_bb-1,j_bb-1,k_bb-1));


% interpolating polynomia contrusction:


B=[ 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    -3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    9,-9,-9, 9, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   -6, 6, 6,-6, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
  -6, 6, 6,-6, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   4,-4,-4, 4, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3,-6,-3, 0, 0, 0, 0, 6,-6, 3,-3, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 3, 3, 0, 0, 0, 0,-4, 4,-2, 2, 0, 0, 0, 0,-2,-2,-1,-1, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 4, 2, 0, 0, 0, 0,-3, 3,-3, 3, 0, 0, 0, 0,-2,-1,-2,-1, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2,-2,-2, 0, 0, 0, 0, 2,-2, 2,-2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0;
  -3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
   9,-9, 0, 0,-9, 9, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
-6, 6, 0, 0, 6,-6, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0, 0, 0,-1, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9,-9, 0, 0,-9, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 3, 0, 0,-6,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6,-6, 0, 0, 3,-3, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3,-3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 4, 0, 0,-2, 2, 0, 0,-2,-2, 0, 0,-1,-1, 0, 0;
 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 9, 0,-9, 0,-9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0,-6, 0,-3, 0, 6, 0,-6, 0, 3, 0,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0;
-27,27,27,-27,27,-27,-27,27,-18,-9,18, 9,18, 9,-18,-9,-18,18,-9, 9,18,-18, 9,-9,-18,18,18,-18,-9, 9, 9,-9,-12,-6,-6,-3,12, 6, 6, 3,-12,-6,12, 6,-6,-3, 6, 3,-12,12,-6, 6,-6, 6,-3, 3,-8,-4,-4,-2,-4,-2,-2,-1;
18,-18,-18,18,-18,18,18,-18, 9, 9,-9,-9,-9,-9, 9, 9,12,-12, 6,-6,-12,12,-6, 6,12,-12,-12,12, 6,-6,-6, 6, 6, 6, 3, 3,-6,-6,-3,-3, 6, 6,-6,-6, 3, 3,-3,-3, 8,-8, 4,-4, 4,-4, 2,-2, 4, 4, 2, 2, 2, 2, 1, 1;
-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 0,-3, 0, 3, 0, 3, 0,-4, 0, 4, 0,-2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-2, 0,-1, 0,-1, 0;
18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6, 9,-9, 9,-9,-9, 9,-9, 9,12,-12,-12,12, 6,-6,-6, 6, 6, 3, 6, 3,-6,-3,-6,-3, 8, 4,-8,-4, 4, 2,-4,-2, 6,-6, 6,-6, 3,-3, 3,-3, 4, 2, 4, 2, 2, 1, 2, 1;
-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-6, 6,-6, 6, 6,-6, 6,-6,-8, 8, 8,-8,-4, 4, 4,-4,-3,-3,-3,-3, 3, 3, 3, 3,-4,-4, 4, 4,-2,-2, 2, 2,-4, 4,-4, 4,-2, 2,-2, 2,-2,-2,-2,-2,-1,-1,-1,-1;
 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
-6, 6, 0, 0, 6,-6, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 4,-4, 0, 0,-4, 4, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-6, 6, 0, 0, 6,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4,-2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0,-2,-1, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4,-4, 0, 0,-4, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0,-2,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0;
-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0,-6, 0, 6, 0, 6, 0,-6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-4, 0,-2, 0, 4, 0, 2, 0,-3, 0, 3, 0,-3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0,-2, 0,-1, 0,-2, 0,-1, 0;
18,-18,-18,18,-18,18,18,-18,12, 6,-12,-6,-12,-6,12, 6,12,-12, 6,-6,-12,12,-6, 6, 9,-9,-9, 9, 9,-9,-9, 9, 8, 4, 4, 2,-8,-4,-4,-2, 6, 3,-6,-3, 6, 3,-6,-3, 6,-6, 3,-3, 6,-6, 3,-3, 4, 2, 2, 1, 4, 2, 2, 1;
-12,12,12,-12,12,-12,-12,12,-6,-6, 6, 6, 6, 6,-6,-6,-8, 8,-4, 4, 8,-8, 4,-4,-6, 6, 6,-6,-6, 6, 6,-6,-4,-4,-2,-2, 4, 4, 2, 2,-3,-3, 3, 3,-3,-3, 3, 3,-4, 4,-2, 2,-4, 4,-2, 2,-2,-2,-1,-1,-2,-2,-1,-1;
 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
 0, 0, 0, 0, 0, 0, 0, 0, 4, 0,-4, 0,-4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0,-2, 0,-2, 0, 2, 0,-2, 0, 2, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0;
-12,12,12,-12,12,-12,-12,12,-8,-4, 8, 4, 8, 4,-8,-4,-6, 6,-6, 6, 6,-6, 6,-6,-6, 6, 6,-6,-6, 6, 6,-6,-4,-2,-4,-2, 4, 2, 4, 2,-4,-2, 4, 2,-4,-2, 4, 2,-3, 3,-3, 3,-3, 3,-3, 3,-2,-1,-2,-1,-2,-1,-2,-1;
 8,-8,-8, 8,-8, 8, 8,-8, 4, 4,-4,-4,-4,-4, 4, 4, 4,-4, 4,-4,-4, 4,-4, 4, 4,-4,-4, 4, 4,-4,-4, 4, 2, 2, 2, 2,-2,-2,-2,-2, 2, 2,-2,-2, 2, 2,-2,-2, 2,-2, 2,-2, 2,-2, 2,-2, 1, 1, 1, 1, 1, 1, 1, 1];


F=[f_000;f_100;f_010;f_110;f_001;f_101;f_011;f_111;fx_000;fx_100;fx_010;fx_110;fx_001;fx_101;fx_011;fx_111;fy_000;fy_100;fy_010;fy_110;fy_001;fy_101;fy_011;fy_111;fz_000;fz_100;fz_010;fz_110;fz_001;fz_101;fz_011;fz_111;fxy_000;fxy_100;fxy_010;fxy_110;fxy_001;fxy_101;fxy_011;fxy_111;fxz_000;fxz_100;fxz_010;fxz_110;fxz_001;fxz_101;fxz_011;fxz_111;fyz_000;fyz_100;fyz_010;fyz_110;fyz_001;fyz_101;fyz_011;fyz_111;fxyz_000;fxyz_100;fxyz_010;fxyz_110;fxyz_001;fxyz_101;fxyz_011;fxyz_111];

C=B*F;

A=zeros(4,4,4);

for i_m=0:3
    for j_m=0:3
        for k_m=0:3
            
            A(i_m+1,j_m+1,k_m+1)=C(1+i_m+4*j_m+16*k_m);
            
            
        end
    end
end 


Phi_int=Phi_3d(i_bottom,j_bottom,k_bottom);

f_x=0;
f_y=0;
f_z=0;

H_xx=0;
H_yy=0;
H_zz=0;
H_xy=0;
H_yx=0;
H_xz=0;
H_zx=0;
H_yz=0;
H_zy=0;

%gradient & hessian computation g= SUM_jSUM_ia_ij*((x-x_0)/step_grid)^i*((y-y_0)/step_grid)^j;;

        for i=0:3
            for j=0:3
                for k=0:3
                    
                Phi_int=Phi_int+A(i+1,j+1,k+1)*((x-x_bottom)/step_grid)^i*((y-y_bottom)/step_grid)^j*((z-z_bottom)/step_grid)^k;
                
                f_x=f_x+A(i+1,j+1,k+1)*(i/step_grid)*((x-x_bottom)/step_grid)^(i-1)*((y-y_bottom)/step_grid)^j*((z-z_bottom)/step_grid)^k;
                f_y=f_y+A(i+1,j+1,k+1)*(j/step_grid)*((x-x_bottom)/step_grid)^i*((y-y_bottom)/step_grid)^(j-1)*((z-z_bottom)/step_grid)^k;
                f_z=f_z+A(i+1,j+1,k+1)*(k/step_grid)*((x-x_bottom)/step_grid)^i*((y-y_bottom)/step_grid)^j*((z-z_bottom)/step_grid)^(k-1);
                
                H_xx=H_xx+A(i+1,j+1,k+1)*(i/step_grid)*((i-1)/step_grid)*((x-x_bottom)/step_grid)^(i-2)*((y-y_bottom)/step_grid)^j*((z-z_bottom)/step_grid)^k;
                H_xy=H_xy+A(i+1,j+1,k+1)*(i/step_grid)*(j/step_grid)*((x-x_bottom)/step_grid)^(i-1)*((y-y_bottom)/step_grid)^(j-1)*((z-z_bottom)/step_grid)^k;
                H_xz=H_xz+A(i+1,j+1,k+1)*(i/step_grid)*(k/step_grid)*((x-x_bottom)/step_grid)^(i-1)*((y-y_bottom)/step_grid)^j*((z-z_bottom)/step_grid)^(k-1);
                H_yx=H_yx+A(i+1,j+1,k+1)*(i/step_grid)*(j/step_grid)*((x-x_bottom)/step_grid)^(i-1)*((y-y_bottom)/step_grid)^(j-1)*((z-z_bottom)/step_grid)^k;
                H_zx=H_zx+A(i+1,j+1,k+1)*(i/step_grid)*(k/step_grid)*((x-x_bottom)/step_grid)^(i-1)*((y-y_bottom)/step_grid)^j*((z-z_bottom)/step_grid)^(k-1);
                H_yy=H_yy+A(i+1,j+1,k+1)*(j/step_grid)*((j-1)/step_grid)*((x-x_bottom)/step_grid)^i*((y-y_bottom)/step_grid)^(j-2)*((z-z_bottom)/step_grid)^k;
                H_yz=H_yz+A(i+1,j+1,k+1)*(k/step_grid)*(j/step_grid)*((x-x_bottom)/step_grid)^i*((y-y_bottom)/step_grid)^(j-1)*((z-z_bottom)/step_grid)^(k-1);
                H_zy=H_zy+A(i+1,j+1,k+1)*(k/step_grid)*(j/step_grid)*((x-x_bottom)/step_grid)^i*((y-y_bottom)/step_grid)^(j-1)*((z-z_bottom)/step_grid)^(k-1);
                H_zz=H_zz+A(i+1,j+1,k+1)*(k/step_grid)*((k-1)/step_grid)*((x-x_bottom)/step_grid)^i*((y-y_bottom)/step_grid)^j*((z-z_bottom)/step_grid)^(k-2);
                end 
                
            end
        end
        
elseif Interp_type==1
    
    
K=Phi_3d(i_bottom,j_bottom,k_bottom);
Ki=Phi_3d(i_top,j_bottom,k_bottom);
Kj=Phi_3d(i_bottom,j_top,k_bottom);
Kk=Phi_3d(i_bottom,j_bottom,k_top);
Kij=Phi_3d(i_top,j_top,k_bottom);
Kik=Phi_3d(i_top,j_bottom,k_top);
Kjk=Phi_3d(i_bottom,j_top,k_top);
Kijk=Phi_3d(i_top,j_top,k_top);  
  
Phi_int=(1/step_grid)^3*(Kijk*(x - x_bottom)*(y - y_bottom)*(z - z_bottom) + Kjk*(-x + x_top)*(y - y_bottom)*(z - z_bottom) + Kik*(x - x_bottom)*(-y + y_top)*(z - z_bottom) + Kk*(-x + x_top)*(-y + y_top)*(z - z_bottom) + Kij*(x - x_bottom)*(y - y_bottom)*(-z + z_top) + Kj*(-x + x_top)*(y - y_bottom)*(-z + z_top) + Ki*(x - x_bottom)*(-y + y_top)*(-z + z_top) + K*(-x + x_top)*(-y + y_top)*(-z + z_top));

f_x= (1/step_grid)^3*(Kijk*(y - y_bottom)*(z - z_bottom) - Kjk*(y - y_bottom)*(z - z_bottom) + Kk*(y - y_top)*(z - z_bottom) + Kik*(-y + y_top)*(z - z_bottom) + Kj*(y - y_bottom)*(z - z_top) + Ki*(y - y_top)*(z - z_top) + Kij*(y - y_bottom)*(-z + z_top) - K*(-y + y_top)*(-z + z_top));
f_y= (1/step_grid)^3*(Kijk*(x - x_bottom)*(z - z_bottom) - Kik*(x - x_bottom)*(z - z_bottom) + Kk*(x - x_top)*(z - z_bottom) + Kjk*(-x + x_top)*(z - z_bottom) + Ki*(x - x_bottom)*(z - z_top) + Kj*(x - x_top)*(z - z_top) + Kij*(x - x_bottom)*(-z + z_top) - K*(-x + x_top)*(-z + z_top));
f_z= (1/step_grid)^3*(Kijk*(x - x_bottom)*(y - y_bottom) - Kij*(x - x_bottom)*(y - y_bottom) + Kj*(x - x_top)*(y - y_bottom) + Kjk*(-x + x_top)*(y - y_bottom) + Ki*(x - x_bottom)*(y - y_top) + Kk*(x - x_top)*(y - y_top) + Kik*(x - x_bottom)*(-y + y_top) - K*(-x + x_top)*(-y + y_top));

H_xx=0;
H_yy=0;
H_zz=0;
H_xy=0;
H_yx=0;
H_xz=0;
H_zx=0;
H_yz=0;
H_zy=0;

    
    

        
elseif Interp_type==0

 %nearest-neighbor interpolation:
   
NN_guess=[norm(r'-[R1_boundedx(i_bottom),R1_boundedy(j_bottom),R1_boundedz(k_bottom)]); %[ 0 0 0]
          norm(r'-[R1_boundedx(i_top),R1_boundedy(j_bottom),R1_boundedz(k_bottom)]);    %[ 1 0 0]
          norm(r'-[R1_boundedx(i_bottom),R1_boundedy(j_top),R1_boundedz(k_bottom)]);    %[ 0 1 0]
          norm(r'-[R1_boundedx(i_bottom),R1_boundedy(j_bottom),R1_boundedz(k_top)]);    %[ 0 0 1]
          norm(r'-[R1_boundedx(i_top),R1_boundedy(j_top),R1_boundedz(k_bottom)]);       %[ 1 1 0]
          norm(r'-[R1_boundedx(i_top),R1_boundedy(j_bottom),R1_boundedz(k_top)]);       %[ 1 0 1]
          norm(r'-[R1_boundedx(i_bottom),R1_boundedy(j_top),R1_boundedz(k_top)]);       %[ 0 1 1]
          norm(r'-[R1_boundedx(i_top),R1_boundedy(j_top),R1_boundedz(k_top)])];         %[ 1 0 1]
    
    
[~,NN]=min(NN_guess);  


if NN==1 %[ 0 0 0]
    
i_bb=i_bottom;
j_bb=j_bottom;
k_bb=k_bottom;

i_tt=i_top;
j_tt=j_top;
k_tt=k_top;

               Phi_int=Phi_3d(i_bottom,j_bottom,k_bottom);
    
elseif NN==2  %[ 1 0 0]
    
    
i_bb=i_top;
j_bb=j_bottom;
k_bb=k_bottom;
    
i_tt=i_top+1;
j_tt=j_top;
k_tt=k_top;
               
                Phi_int=Phi_3d(i_top,j_bottom,k_bottom);

elseif NN==3 %[ 0 1 0]
    
i_bb=i_bottom;
j_bb=j_top; 
k_bb=k_bottom;
    
i_tt=i_top;
j_tt=j_top+1;    
k_tt=k_top;
               
               Phi_int=Phi_3d(i_bottom,j_top,k_bottom);
               
elseif NN==4  %[ 0 0 1]
    
    
i_bb=i_bottom;
j_bb=j_bottom;
k_bb=k_top;

i_tt=i_top;
j_tt=j_top;
k_tt=k_top+1;
 
               Phi_int=Phi_3d(i_bottom,j_bottom,k_top);
elseif NN==5  %[ 1 1 0]
    
i_bb=i_top;
j_bb=j_top;
k_bb=k_bottom;
    
i_tt=i_top+1;
j_tt=j_top+1;
k_tt=k_top;
            Phi_int=Phi_3d(i_top,j_top,k_bottom);
elseif NN==6  %[ 1 0 1]

i_bb=i_top;
j_bb=j_bottom;
k_bb=k_top;
    
i_tt=i_top+1;
j_tt=j_top;
k_tt=k_top+1;
             Phi_int=Phi_3d(i_top,j_bottom,k_top);

elseif NN==7 %[0 1 1]
    
i_bb=i_bottom;
j_bb=j_top;
k_bb=k_top;

i_tt=i_top;
j_tt=j_top+1;
k_tt=k_top+1;    

               Phi_int=Phi_3d(i_bottom,j_top,k_top);

elseif NN==8 %[1 1 1]

i_bb=i_top;
j_bb=j_top;
k_bb=k_top;

i_tt=i_top+1;
j_tt=j_top+1;
k_tt=k_top+1;    

            Phi_int=Phi_3d(i_top,j_top,k_top);

end 

            

               
               

f_x=1/(2*step_grid)*(Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb-1,j_bb,k_bb));
f_y=1/(2*step_grid)*(Phi_3d(i_bb,j_tt,k_bb)-Phi_3d(i_bb,j_bb-1,k_bb));
f_z=1/(2*step_grid)*(Phi_3d(i_bb,j_bb,k_tt)-Phi_3d(i_bb,j_bb,k_bb-1));

H_xx=1/(step_grid^2)*( Phi_3d(i_tt,j_bb,k_bb)-2*Phi_3d(i_bb,j_bb,k_bb)+Phi_3d(i_bb-1,j_bb,k_bb));
H_xy=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb-1,k_bb)- Phi_3d(i_bb-1,j_bb,k_bb) - Phi_3d(i_bb,j_bb-1,k_bb)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)+ Phi_3d(i_tt,j_tt,k_bb));
H_yx=H_xy;
H_xz=1/(2*step_grid^2)*( Phi_3d(i_bb-1,j_bb,k_bb-1)- Phi_3d(i_bb-1,j_bb,k_bb) - Phi_3d(i_bb,j_bb,k_bb-1)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_tt,j_bb,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_tt,j_bb,k_tt));
H_zx=H_xz;
H_yy=1/(step_grid^2)*( Phi_3d(i_bb,j_tt,k_bb)-2*Phi_3d(i_bb,j_bb,k_bb)+Phi_3d(i_bb,j_bb-1,k_bb));
H_yz=1/(2*step_grid^2)*( Phi_3d(i_bb,j_bb-1,k_bb-1)- Phi_3d(i_bb,j_bb-1,k_bb) - Phi_3d(i_bb,j_bb,k_bb-1)+2*Phi_3d(i_bb,j_bb,k_bb)-Phi_3d(i_bb,j_tt,k_bb)-Phi_3d(i_bb,j_bb,k_tt)+ Phi_3d(i_bb,j_tt,k_tt));
H_zy=H_yz;
H_zz=1/(step_grid^2)*( Phi_3d(i_bb,j_bb,k_tt)-2*Phi_3d(i_bb,j_bb,k_bb)+Phi_3d(i_bb,j_bb,k_bb-1));

    
    
    
end

  
f=[f_x;f_y;f_z];

H=[H_xx, H_xy,  H_xz;
   H_yx, H_yy,  H_yz;
   H_zx, H_zy,  H_zz];



return