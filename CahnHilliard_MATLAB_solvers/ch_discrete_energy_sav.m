function E = ch_discrete_energy_sav(phi,hxy,eps2,gamma0)
    [gridx,gridy] = size(phi);
    a = hxy*sum(sum(f_SAV(phi,0))); %Calculate chemical free energy; gamma0 is set to 0
    sum_i = 0; % Initialize interfacial free energy in x
    for i = 1:gridx-1
        for j = 1:gridy
            sum_i = sum_i + (phi(i+1,j)-phi(i,j))^2;
        end
    end
    sum_j = 0; % Initialize interfacial free energy in y
    for i = 1:gridx
        for j = 1:gridy-1
            sum_j = sum_j + (phi(i,j+1)-phi(i,j))^2;
        end
    end
    E = a + 0.5*eps2*(sum_i+sum_j);
end

% Unused code from Min-Jhe
% function E = ch_discrete_energy_sav(phi,hxy,eps2,k2,gamma0,r,C0)
%     % % % E1 = hxy*fft2(f(phi,0)); %Calculate chemical free energy
%     % % E1 = hxy*fft2(f(phi,gamma0)); %Calculate chemical free energy
%     % % E1 = E1(1,1);

%     % E = 0.5 * fft2(phi .*(-eps2*real(ifft2(k2 .* fft2(phi)))));
%     % E = r^2 + E(1,1);

%     E_gamma = hxy*fft2(gamma0/2*phi.^2);
%     E_gamma = E_gamma(1,1);

%     [gridx,gridy] = size(phi);
%     a = r^2-C0-(gamma0^2+2*gamma0)/4; 
%     sum_i = 0; % Initialize interfacial free energy in x
%     for i = 1:gridx-1
%         for j = 1:gridy
%             sum_i = sum_i + (phi(i+1,j)-phi(i,j))^2;
%         end
%     end
%     sum_j = 0; % Initialize interfacial free energy in y
%     for i = 1:gridx
%         for j = 1:gridy-1
%             sum_j = sum_j + (phi(i,j+1)-phi(i,j))^2;
%     end
%     E = a + 0.5*eps2*(sum_i+sum_j)+E_gamma;

% end

% function E = ch_discrete_energy(phi, hxy, eps2, gamma0)
%     % Get grid dimensions
%     [gridx, gridy] = size(phi);
    
%     % Calculate chemical free energy (unchanged)
%     a = hxy * sum(sum(f(phi, 0)));
    
%     % Compute discrete Laplacian with periodic boundaries
%     phi_right = circshift(phi, [0, 1]);   % Shift right (y-direction)
%     phi_left = circshift(phi, [0, -1]);   % Shift left (y-direction)
%     phi_up = circshift(phi, [-1, 0]);     % Shift up (x-direction)
%     phi_down = circshift(phi, [1, 0]);    % Shift down (x-direction)
%     Lap_phi = phi_right + phi_left + phi_up + phi_down - 4 * phi;
%     % Note: Assuming h=1; if h≠1, divide Lap_phi by h^2 and adjust hxy accordingly
    
%     % Compute interfacial free energy using integration by parts
%     interfacial_energy = - (eps2 / 2) * hxy * sum(sum(phi .* Lap_phi));
    
%     % Total energy
%     E = a + interfacial_energy;
% end