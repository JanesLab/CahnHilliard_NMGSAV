function lap_a = ch_laplace(a,nx,ny,xright,xleft,yright,yleft,boundary)
%This function calculates the discrete 2D Laplacian for matrix a.
%
%INPUTS
    %a = Matrix to be operated on.
    %nx = Number of grid points in x.
    %ny = Number of grid points in y.
    %xright = Rightmost grid point in x.
    %xleft = Leftmost grid point in x.
    %yright = Rightmost gridpoint in y.
    %yleft = Leftmost gridpoint in y.
    %boundary = 'periodic' or 'neumann' boundary conditions.
%
%OUTPUT
    %lap_a = discrete 2D Laplacian of a.

lap_a = zeros(nx,ny); %Initialize Laplacian
ht2 = ((xright-xleft)/nx)*((yright-yleft)/ny); %Define grid size
%Determine boundary condition
periodic = strcmpi(boundary,'periodic');
neumann = strcmpi(boundary,'neumann');

for i = 1:nx
    for j = 1:ny
        if i > 1
            dadx_L = a(i,j)-a(i-1,j);
        elseif neumann
            dadx_L = 0;
        elseif periodic
            dadx_L = a(i,j)-a(nx,j); %changed
        end
        if i < nx
            dadx_R = a(i+1,j)-a(i,j);
        elseif neumann
            dadx_R = 0;
        elseif periodic
            dadx_R = a(1,j)-a(i,j);
        end
        if j > 1
            dady_B = a(i,j)-a(i,j-1);
        elseif neumann
            dady_B = 0;
        elseif periodic
            dady_B = a(i,j)-a(i,ny); %changed
        end
        if j < ny
            dady_T = a(i,j+1)-a(i,j);
        elseif neumann
            dady_T = 0;
        elseif periodic
            dady_T = a(i,1)-a(i,j);
        end
        lap_a(i,j) = (dadx_R-dadx_L + dady_T-dady_B)/ht2;
    end
end

end