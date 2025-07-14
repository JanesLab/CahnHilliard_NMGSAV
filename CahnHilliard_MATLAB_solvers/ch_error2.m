function res2 = ch_error2(phi_oldf,phi_newf,muf,nxf,nyf,xr,xl,yr,yl,dtf,bc)
%This function computes the 2D residual for phi.
%
%INPUTS
    %phi_oldf = Prior chemical state.
    %phi_newf = Next chemical state.
    %muf = Chemical potential.
    %nxf = Number of grid points in x.
    %nyf = Number of grid points in y.
    %xr = Rightmost grid point in x.
    %xl = Leftmost grid point in x.
    %yr = Rightmost gridpoint in y.
    %yl = Leftmost gridpoint in y.
    %dtf = Time step.
    %bc = 'periodic' or 'neumann' boundary conditions.
%
%OUTPUT
    %res2 = Normalized residual

rr = muf - phi_oldf; %Calculate the starting residual
sor = ch_laplace(rr,nxf,nyf,xr,xl,yr,yl,bc); %Calculate the source term from rr
rr = sor - (phi_newf-phi_oldf)/dtf; %Update the residual
res2 = sqrt(sum(sum(rr.^2))/(nxf*nyf)); %Take Frobenius norm of sum of squared error

end
