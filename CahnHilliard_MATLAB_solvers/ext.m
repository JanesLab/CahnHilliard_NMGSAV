% Local function for "flip" extension
function x_ext = ext(x)
    % Mirroring extension: takes x(nx, ny) -> x_ext(2*nx, 2*ny)
        [nx, ny] = size(x);
        x_ext = zeros(2*nx, 2*ny);
    
        % Original block
        x_ext(1:nx, 1:ny) = x;
    
        % Flip horizontally
        x_ext(1:nx, ny+1:2*ny) = x(:, end:-1:1);
    
        % Flip vertically
        x_ext(nx+1:2*nx, 1:ny) = x(end:-1:1, :);
    
        % Flip both
        x_ext(nx+1:2*nx, ny+1:2*ny) = x(end:-1:1, end:-1:1);
    end