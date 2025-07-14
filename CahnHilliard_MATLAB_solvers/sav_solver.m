function [phi_new,r_new] = sav_solver(phi_old, phi_prev, r_old, ...
    hx,hy,k2,k4,dt,epsilon2, ...
    boundary,C0,gamma0,eta,xi_flag,i)
%This function uses the sav method to solve the 
%Cahn-Hilliard equation for the next time step.
%
%INPUTS
    % phi_old  = Prior chemical state.
    % phi_prev = Previous chemical state.
    % r_old    = Prior SAV parameter.
    % hx, hy   = Grid spacings in x and y directions.
    % k2, k4   = Laplacian and biharmonic operators.
    % dt       = Time step.
    % epsilon2 = Square of the interface width parameter.
    % boundary = 'periodic' or 'neumann'.
    % C0       = Regularization parameter.
    % gamma0   = Stabilization parameter.
    % eta      = Relaxation parameter.
    % xi_flag  = Optional flag (0 or 1); if 1, xi is set to 1 (no relaxation).
    % i        = Current time step index.
%
%OUTPUT
    %phi_new = Next chemical state.
    %r_new   = Next sav state.
        
    %% CN(Crank-Nicolson) time-stepping
    phi0 = phi_old;
    r0   = r_old;

    phi0_df   = df(phi0,gamma0);        % df at phi0
    Lap_dfphi0 = Lap_SAV(phi0_df, k2, boundary);                    % Lap of df(phi0)
    if i == 1
        phi_bar = A_inv_CN(phi0 + dt/2 * Lap_dfphi0, dt, k2, k4, gamma0, epsilon2, boundary);
    elseif i>=2
        % phi_bar = A_inv_CN(phi0 + dt/2 * Lap_dfphi0, dt, k2, k4, gamma0, epsilon2,  boundary);
        phi_bar = 1.5*phi_old - 0.5*phi_prev;
        phi_bar = max(-1, min(1, phi_bar)); 
        % fprintf('Step %d:  max|phi_bar|=%.3f,  diff_prev=%.3e\n', ...
        %          i, max(abs(phi_bar(:))), norm(phi_old(:)-phi_prev(:)));
    end



    % Step 1
    b = b_fun(phi_bar,hx,hy,C0,gamma0);
    g = g_fun_CN(phi0,r0,b,dt,hx,hy,epsilon2,gamma0,C0,k2,boundary);

    AiLb = A_inv_CN(Lap_SAV(b,k2,boundary),dt,k2,k4,gamma0,epsilon2,boundary);
    Aig  = A_inv_CN(g,dt,k2,k4,gamma0,epsilon2,boundary);

    gamma = -fft2(b.*AiLb);
    gamma = gamma(1,1)*hx*hy;

    % Step 2
    bphi = fft2(b.*Aig);
    bphi = bphi(1,1)*hx*hy/(1+dt/4*gamma);

    % Step 3
    phi_new = dt/4*bphi.*AiLb + Aig;
    r_new   = r_fun(phi_new, phi_old, r0, b, hx, hy, C0, dt, gamma0);

    % Calculate a, b, c
        % Q_phi_new
            E1_new = fft2(f_SAV(phi_new, gamma0));
            E1_new = E1_new(1,1)*hx*hy;
            Q_phi_new = sqrt(E1_new+C0);
        % q_tilde
            q_tilde = r_new;
        % mu_half and Lap_mu_half
            phi_half = (phi_new + phi_old) / 2;
            r_half   = (r_new + r_old) / 2; 
            mu_half = Lap_SAV(phi_half, k2, boundary) + r_half .* b;
            Lap_mu_half = Lap_SAV(mu_half, k2, boundary);
        % muGmu_half
            muGmu_half = fft2(mu_half .* (-Lap_mu_half));
            muGmu_half = muGmu_half(1,1)*hx*hy;
        % a, b, c
            a = (q_tilde - Q_phi_new)^2;
            b = 2 * (q_tilde - Q_phi_new) * Q_phi_new;
            c = Q_phi_new^2 - q_tilde^2 - dt * eta * muGmu_half;

    % Calculate relaxation parameter xi based on xi_flag
        if xi_flag == 0
            xi = 1; % No relaxation when xi_flag is 1
        elseif xi_flag == 1
            if a > 0
                discriminant = b^2 - 4 * a * c; % Calculate discriminant
                if discriminant >= 0
                    xi = (-b - sqrt(discriminant)) / (2 * a);
                    xi = max(0, min(xi, 1)); % Restrict xi to the range [0,1]
                else
                    xi = 1; % When discriminant < 0, choose xi = 1
                end
            else
                % Special case when a = 0
                if b ~= 0
                    xi = -c / b;
                    xi = max(0, min(xi, 1)); % Restrict xi to the range [0,1]
                else
                    if c <= 0
                        xi = 1; % When c <= 0, choose xi = 1
                    else
                        xi = 0; % When c > 0, no solution, choose xi = 0
                    end
                end
            end
        end
    
    % Update r_new
        r_new = xi * q_tilde + (1 - xi) * Q_phi_new;
end