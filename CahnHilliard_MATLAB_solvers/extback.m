
% Local function for "flip" extension back
function x_back = extback(x_ext)
    % Shrinks from 2*nx x 2*ny back to nx x ny (upper-left block)
        [nx_ext, ny_ext] = size(x_ext);
        nx = nx_ext/2;
        ny = ny_ext/2;
        x_back = x_ext(1:nx, 1:ny);
end