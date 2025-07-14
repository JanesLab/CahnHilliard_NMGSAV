function Lphi = Lap_SAV(phi,k2,boundary)
% This function computes the Laplacian of phi using the Fourier transform.
% Inputs:
%   phi      - Current chemical state.
%   k2       - Laplacian operator.
%   boundary  - Boundary condition ('periodic' or 'neumann').
% Outputs:
%   Lphi     - Laplacian of phi.

% Check boundary condition
if strcmp(boundary, 'periodic')
    % Apply periodic boundary conditions
    Lphi = real(ifft2(k2 .* fft2(phi)));
elseif strcmp(boundary, 'neumann')
    % Apply Neumann boundary conditions
    % Use a filtered version of phi to apply cosine transform
    Lphi = real(ifft2(k2 .* fft2_filtered(phi)));
end