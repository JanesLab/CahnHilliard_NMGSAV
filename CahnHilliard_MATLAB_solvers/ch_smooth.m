function [phi_new,mu_new] = ch_smooth(phi_old,mu_old,s_phi,s_mu,c_relax,...
    nx,ny,xright,xleft,yright,yleft,dt,epsilon2,boundary)
%This smooths chemical state and chemical potential maps by the number of 
%cycles specified by c_relax.
%
%INPUTS
    %phi_old = Initial chemical state.
    %mu_old = Initial chemical potential.
    %s_phi = Source term for chemical state (see ζ term in Lee et al., Mathematics 2020).
    %s_mu = Source term for chemical postential (see Ψ term in Lee et al., Mathematics 2020).
    %c_relax = Number of smoothing relaxations.
    %nx = Number of grid points in x.
    %ny = Number of grid points in y.
    %xright = Rightmost grid point in x.
    %xleft = Leftmost grid point in x.
    %yright = Rightmost gridpoint in y.
    %yleft = Leftmost gridpoint in y.
    %dt = Time step.
    %epsilon2 = Squared interfacial transition distance (see ϵ^2 term in Lee et al., Mathematics 2020).
    %boundary = 'periodic' or 'neumann' boundary conditions.
%
%OUTPUT
    %phi_new = smoothed chemical state
    %mu_new = smoothed chemical potential

ht2 = ((xright-xleft)/nx)*((yright-yleft)/ny); %Define grid size
a = zeros(2); %Initialize Jacobian
f = zeros(2,1); %Initialize source function
%Determine boundary condition
periodic = strcmpi(boundary,'periodic');
neumann = strcmpi(boundary,'neumann');

phi_new = phi_old;
mu_new = mu_old;

for iter = 1:c_relax
    for i = 1:nx
        for j = 1:ny
            if neumann
                if i > 1 && i < nx
                    x_fac = 2;
                else
                    x_fac = 1;
                end
                if j > 1 && j < ny
                    y_fac = 2;
                else
                    y_fac = 1;
                end
            end
            a(1,1) = 1/dt;
            if neumann
                a(1,2) = (x_fac+y_fac)/ht2;
                a(2,1) = -epsilon2*(x_fac+y_fac)/ht2 - 3*phi_new(i,j)^2;
            elseif periodic
                a(1,2) = 4/ht2;
                a(2,1) = -4*epsilon2/ht2 - 3*phi_new(i,j)^2;
            end
            a(2,2) = 1;
            f(1) = s_phi(i,j);
            f(2) = s_mu(i,j) - 2*phi_new(i,j)^3;
            if i > 1
                f(1) = f(1) + mu_new(i-1,j)/ht2;
                f(2) = f(2) - epsilon2*phi_new(i-1,j)/ht2;
            elseif periodic
                f(1) = f(1) + mu_new(nx,j)/ht2; %changed
                f(2) = f(2) - epsilon2*phi_new(nx,j)/ht2; %changed
            end
            if i < nx
                f(1) = f(1) + mu_new(i+1,j)/ht2;
                f(2) = f(2) - epsilon2*phi_new(i+1,j)/ht2;
            elseif periodic
                f(1) = f(1) + mu_new(1,j)/ht2; %changed
                f(2) = f(2) - epsilon2*phi_new(1,j)/ht2; %changed
            end
            if j > 1
                f(1) = f(1) + mu_new(i,j-1)/ht2;
                f(2) = f(2) - epsilon2*phi_new(i,j-1)/ht2;
            elseif periodic
                f(1) = f(1) + mu_new(i,ny)/ht2; %changed
                f(2) = f(2) - epsilon2*phi_new(i,ny)/ht2;%changed
            end
            if j < ny
                f(1) = f(1) + mu_new(i,j+1)/ht2;
                f(2) = f(2) - epsilon2*phi_new(i,j+1)/ht2;
            elseif periodic
                f(1) = f(1) + mu_new(i,1)/ht2; %changed
                f(2) = f(2) - epsilon2*phi_new(i,1)/ht2; %changed
            end
            phi_new(i,j) = (a(2,2)*f(1)-a(1,2)*f(2))/det(a); %Solve for next phi
            mu_new(i,j) = (-a(2,1)*f(1) + a(1,1)*f(2))/det(a); %Solve for next mu
        end
    end
end

end