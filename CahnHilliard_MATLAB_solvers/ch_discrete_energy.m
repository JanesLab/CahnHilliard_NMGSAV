function E = ch_discrete_energy(phi,hxy,eps2)
    [gridx,gridy] = size(phi);
    a = hxy*sum(sum(f(phi))); %Calculate chemical free energy; gamma0 is set to 0
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
