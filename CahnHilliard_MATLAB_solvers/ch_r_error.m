function D = ch_r_error(r,phi,hxy,C0,gamma0)
    % D=(r-sqrt(hxy*sum(sum(f(phi)))))/sqrt(hxy*sum(sum(f(phi))));
    E1 = fft2(f_SAV(phi,gamma0));
    D  = r-sqrt(hxy*E1(1,1)+C0);
    D  = D/sqrt(hxy*E1(1,1)+C0); % Normalize the error
end