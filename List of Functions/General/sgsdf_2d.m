function h=sgsdf_2d(x,y,nx,ny,d,flag_coupling)
% https://www.mathworks.com/matlabcentral/fileexchange/8991-2-d-savitzky-golay-smoothing-and-differentiation-filter
% 2-D Savitzky-Golay Smoothing and Differentiation Filter
% The filter coefficients are calculated from a matrix approach. 
%   h=sgsdf_2d(x,y,nx,ny,d,flag_coupling)
%       x    = x data point, e.g., -3:3 (both odd-length and even-lenght are accepted, symmetry is prefered)
%       y    = y data point, e.g., -2:2  (both odd-length and even-lenght are accepted)
%       nx    =x polynomial order       default=1            
%       ny    =y polynomial order       default=1
%      d    = differentiation order (0<=d<=nx)    default=0=smoothing
%       flag_coupling  = with or without the consideration of the coupling
%                                   terms, between x and y. default=0
%       h   = filter coefficients obtained.
% Example:
%     sgsdf_2d(-2:2,-3:3,1,1,0,0)
%     sgsdf_2d(-2:2,-3:3,2,2,0,1)
%     sgsdf_2d((-3:2)+1/2,-3:3,2,2,1,0)
%     sgsdf_2d((-3:2)+1/2,(-4:3)+1/2,2,2,1,1)
%     sgsdf_2d(-3:3,-4:4,2,3,2,1)
% Author:
%     Jianwen Luo <luojw@ieee.org>
%     11/10/2005
% References:
% [1]	J. W. Luo, K. Ying, P. He, and J. Bai, "Properties of Savitzky-Golay Digital Differentiators," 
% Digital Signal Processing, vol. 15, pp. 122-136, 2005.
% [2]	A. Savitzky and M. J. E. Golay, "Smoothing and Differentiation of Data by Simplified Least Squares Procedures," 
% Analytical Chemistry, vol. 36, pp. 1627-1639, 1964.
% [3]	K. L. Ratzlaff and J. T. Johnson, "Computation of Two-Dimensional Polynomial Least-Squares Convolution Smoothing Integers," 
% Analytical Chemistry, vol. 61, pp. 1303-1305, 1989.
% [4]	J. E. Kuo, H. Wang, and S. Pickup, "Multidimensional Least-Squares Smoothing Using Orthogonal Polynomials," 
% Analytical Chemistry, vol. 63, pp. 630-635, 1991.
% [5] http://research.microsoft.com/users/jckrumm/SavGol/SavGol.htm
% [6]	T. A. Sakharuk, "Computation of Weighting Functions for Smoothing 2-Dimensional Data by Local Polynomial-Approximation Techniques," 
% Analytica Chimica Acta, vol. 249, pp. 331-336, 1991.
% [7]	P. Meer and I. Weiss, "Smoothed Differentiation Filters for Images," 
% Journal of Visual Communication and Image Representation, vol. 3, pp. 58-72, 1992.
% [8]	P. Nikitas and A. Pappa-Louisi, "Comments on the two-dimensional smoothing of data," 
% Analytica Chimica Acta, vol. 415, pp. 117-125, 2000.
if nargin<6
    flag_coupling=0;
end
if nargin<5
    d=0;
end
if nargin<4
    ny=1;
end
if nargin<3
    nx=1;
end
[x,y]=meshgrid(x,y);
[ly,lx]=size(x);
if nx>lx-1 || ny>ly-1
    error('polynomial order too large!');
end
if d>nx
    error('differentiation order too large!');    
end
x=x(:);
y=y(:);
A=[];
if flag_coupling
    for kx=nx:-1:0 %nx,nx-1,...,x,1 (nx+1)
        for ky=ny:-1:0 %ny,ny-1,...,y,1 (ny+1)
            A=[A x.^kx.*y.^ky];
        end
    end     
else
    for k=nx:-1:0
        A=[A x.^k];
    end
    for k=ny:-1:1
        A=[A y.^k];
    end  
end
% A'*A
% inv(A'*A)
h=inv(A'*A)*A';
if flag_coupling
    h=h((nx+1-d)*(ny+1),:);% ?*x 
else    
    h=h(nx+1-d,:);% ?*x
end
h=reshape(h,ly,lx);
h=prod(1:d)*h;
