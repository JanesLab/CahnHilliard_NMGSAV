function E = ch_modified_energy(phi,hxy,eps2,k2,gamma0)
    % E1 = hxy*fft2(f(phi,gamma0)); %Calculate chemical free energy
    E1 = hxy*fft2(f_SAV(phi,0)); %Calculate chemical free energy
    E1 = E1(1,1);

    E = 0.5*eps2*fft2(phi.*(-real(ifft2(k2 .* fft2(phi)))));
    E = E1 + E(1,1);
end