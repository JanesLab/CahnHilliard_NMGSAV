function r0 = r0_fun(phi0,hx,hy,C0,gamma0)
    % fphi = f(phi0,gamma0);
    % r0 = sqrt(hx*hy*sum(sum(fphi)) + C0);
    E1 = fft2(f_SAV(phi0,gamma0));

    r0 = sqrt(hx*hy*E1(1,1) + C0);
end