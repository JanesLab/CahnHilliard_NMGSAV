function [uf_new,wf_new] = nmg_vcycle(uf_new,wf_new,su,sw,...
    nx,ny,xright,xleft,yright,yleft,ilevel,c_relax,dt,epsilon2,n_level, ...
    boundary)
%This function implements the V-cycle algorithm recursively.
%
%INPUTS
    %uf_new = Original chemical state.
    %wf_new = Original chemical potential.
    %su = Source term for chemical state (see ζ term in Lee et al., Mathematics 2020).
    %sw = Source term for chemical postential (see Ψ term in Lee et al., Mathematics 2020).
    %nx = Number of grid points in x.
    %ny = Number of grid points in y.
    %xright = Rightmost grid point in x.
    %xleft = Leftmost grid point in x.
    %yright = Rightmost gridpoint in y.
    %yleft = Leftmost gridpoint in y.
    %ilevel = Current V-cycle level
    %c_relax = Number of smoothing relaxations.
    %dt = Time step.
    %epsilon2 = Squared interfacial transition distance (see ϵ^2 term in Lee et al., Mathematics 2020).
    %n_level = Total number of V-cycles.
    %boundary = 'periodic' or 'neumann' boundary conditions.
%
%OUTPUT
    %uf_new = Updated chemical state.
    %wf_new = Updated chemical potential.

[uf_new,wf_new] = ch_smooth(uf_new,wf_new,su,sw,c_relax,nx,ny, ...
    xright,xleft,yright,yleft,dt,epsilon2,boundary); %Smooth the input data
if ilevel < n_level
    uc_new = restrict(uf_new); %Restrict chemical state
    wc_new = restrict(wf_new); %Restrict chemical potential
    [duc,dwc] = defect(uf_new,wf_new,su,sw,nx,ny,uc_new,wc_new, ...
        floor(nx/2),floor(ny/2),xright,xleft,yright,yleft,dt,epsilon2, ...
        boundary); %Compute defect
    uc_def = uc_new;
    wc_def = wc_new;
    [uc_def,wc_def] = nmg_vcycle(uc_def,wc_def,duc,dwc, ...
        floor(nx/2),floor(ny/2),xright,xleft,yright,yleft,ilevel+1, ...
        c_relax,dt,epsilon2,n_level,boundary); %Recursively call V-cycle
    uc_def = uc_def-uc_new;
    wc_def = wc_def-wc_new;
    uf_def = expand(uc_def); %Expand chemical state
    wf_def = expand(wc_def); %Expand chemical potential
    uf_new = uf_new+uf_def;
    wf_new = wf_new+wf_def;
    [uf_new,wf_new] = ch_smooth(uf_new,wf_new,su,sw,c_relax,nx,ny, ...
        xright,xleft,yright,yleft,dt,epsilon2,boundary); %Post-smoothing step
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mc = restrict(mf)
%This function restricts a matrix twofold in each direction.
%
%INPUTS
    %mf = Original matrix.
%
%OUTPUT
    %mc = Compressed matrix.

[nx,ny] = size(mf);
nxc = floor(nx/2);
nyc = floor(ny/2);
mc = zeros(nxc,nyc);

for i = 1:nxc
    for j = 1:nyc
        mc(i,j) = 0.25 * (mf(2*i,2*j) + mf(2*i-1,2*j) + ...
            mf(2*i,2*j-1) + mf(2*i-1,2*j-1));
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ducf,dwcf] = defect(uf_newf,wf_newf,suf,swf,nxf,nyf, ...
    uc_newf,wc_newf,nxc,nyc,xr,xl,yr,yl,dtf,epsilon2f,bc)
%This function computes the defect within the V-cycle algorithm.
%
%INPUTS
    %uf_newf = Original chemical state.
    %wf_newf = Original chemical potential.
    %suf = Original source term for chemical state (see ζ term in Lee et al., Mathematics 2020)
    %swf = Original source term for chemical postential (see Ψ term in Lee et al., Mathematics 2020)
    %nxf = Original number of grid points in x.
    %nyf = Original number of grid points in y.
    %uc_newf = Compressed chemical state.
    %wc_newf = Compressed chemical potential.
    %nxc = Compressed number of grid points in x.
    %nyc = Compressed number of grid points in y.
    %xr = Rightmost grid point in x.
    %xl = Leftmost grid point in x.
    %yr = Rightmost gridpoint in y.
    %yl = Leftmost gridpoint in y.
    %dtf = Time step.
    %epsilon2f = Squared interfacial transition distance (see ϵ^2 term in Lee et al., Mathematics 2020).
    %bc = 'periodic' or 'neumann' boundary conditions.
%
%OUTPUT
    %ducf = Defect in compressed chemical state.
    %dwcf = Defect in compressed chemical potential.

ruf = uf_newf/dtf - ch_laplace(wf_newf,nxf,nyf,xr,xl,yr,yl,bc);
rwf = wf_newf/dtf - uf_newf.^3 + epsilon2f*ch_laplace(uf_newf,nxf,nyf,xr,xl,yr,yl,bc);
ruc = uc_newf/dtf - ch_laplace(wc_newf,nxc,nyc,xr,xl,yr,yl,bc);
rwc = wc_newf/dtf - uc_newf.^3 + epsilon2f*ch_laplace(uc_newf,nxc,nyc,xr,xl,yr,yl,bc);
ruf = suf-ruf;
rwf = swf-rwf;
ducf = ruc + restrict(ruf);
dwcf = rwc + restrict(rwf);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mf = expand(mc)
%This function extends a matrix twofold in each direction.
%
%INPUTS
    %mc = Original matrix.
%
%OUTPUT
    %mf = Extended matrix.

[nx,ny] = size(mc);
nxe = 2*nx;
nye = 2*ny;
mf = zeros(nxe,nye);

for i = 1:nx
    for j = 1:ny
        mf(2*i,2*j) = mc(i,j);
        mf(2*i-1,2*j) = mc(i,j);
        mf(2*i,2*j-1) = mc(i,j);
        mf(2*i-1,2*j-1) = mc(i,j);
    end
end

end