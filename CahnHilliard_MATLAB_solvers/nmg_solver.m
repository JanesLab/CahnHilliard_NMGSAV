function phi_new = nmg_solver(phi_old,phi_new,mu,nx,ny, ...
    xright,xleft,yright,yleft,c_relax,dt,epsilon2,n_level,max_iter,tol, ...
    boundary,printres, pathname)
%This function uses the nonlinear multigrid method to solve the 
%Cahn-Hilliard equation for the next time step.
%
%INPUTS
    %phi_old = Prior chemical state.
    %phi_new = Next chemical state.
    %mu = Chemical potential
    %nx = Number of grid points in x.
    %ny = Number of grid points in y.
    %xright = Rightmost grid point in x.
    %xleft = Leftmost grid point in x.
    %yright = Rightmost gridpoint in y.
    %yleft = Leftmost gridpoint in y.
    %n_level = Number V-cycle levels.
    %c_relax = Number of smoothing relaxations.
    %dt = Time step.
    %epsilon2 = Squared interfacial transition distance (see ϵ^2 term in Lee et al., Mathematics 2020)
    %n_level = Total number of V-cycles.
    %max_iter = Maximum number of iterations for numerical convergence
    %tol = Tolerance for numerical convergence
    %boundary = 'periodic' or 'neumann' boundary conditions.
    %printres = Logical to print residuals to file.
    %pathname = Name of the path to which residuals are printed. Default = 'cd'. May include a prefix for the filename.
%
%OUTPUT
    %phi_new = Next chemical state.

iter = 0;
resid2 = 1;

%Calculate source terms
s_mu = zeros(nx,ny);
phi_lap = ch_laplace(phi_old,nx,ny,xright,xleft,yright,yleft,boundary);
s_phi = phi_old/dt - phi_lap;

while resid2 > tol && iter < max_iter
    [phi_new,mu] = nmg_vcycle(phi_new,mu,s_phi,s_mu,nx,ny, ...
        xright,xleft,yright,yleft,1,c_relax,dt,epsilon2,n_level,boundary); %V-cycle starting a level 1
    resid2 = ch_error2(phi_old,phi_new,mu,nx,ny,xright,xleft,yright,yleft,dt,boundary);
    if printres
        writematrix(resid2,strcat(pathname,'_residuals.csv'),'WriteMode','append');
    end
    iter = iter+1;
end

end