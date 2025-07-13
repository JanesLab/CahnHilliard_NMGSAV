function Ai = A_inv_CN(phi,dt,k2,k4,gamma0,epsilon2,boundary)
    % A_inv_CN computes the inverse of the operator A in the Cahn-Hilliard equation
    % using the Crank-Nicolson method.
    % Inputs:
    %   phi      - Current chemical state.
    %   dt       - Time step.
    %   k2       - Laplacian operator.
    %   k4       - Biharmonic operator.
    %   gamma0   - Stabilization parameter.
    %   epsilon2 - Square of the interface width parameter.
    %   boundary  - Boundary condition ('periodic' or 'neumann').
    % Outputs:
    %   Ai       - Inverse operator applied to phi.

    % Check boundary condition
    if strcmp(boundary, 'periodic')
        % Apply periodic boundary conditions
        denom = 1 + (dt/2*epsilon2*k4 - dt/2*gamma0*k2);
        Ai = real(ifft2(fft2(phi)./denom));
    elseif strcmp(boundary, 'neumann')
        % Apply Neumann boundary conditions
        denom = 1 + (dt/2*epsilon2*k4 - dt/2*gamma0*k2);
        Ai = real(ifft2(fft2_filtered(phi)./denom));
    else
        error('Invalid boundary condition. Use ''periodic'' or ''neumann''.');
    end

end