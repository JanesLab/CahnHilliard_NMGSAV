function r = r_fun(phi,phi0,r0,b,hx,hy,C0,dt,gamma0)
    bphi0 = fft2(b.*phi0);
    bphi0 = hx*hy*bphi0(1,1);
    bphi  = fft2(b.*phi);
    bphi  = hx*hy*bphi(1,1);

    E1 = fft2(f_SAV(phi0,gamma0));
    r = r0 + 1/2*(bphi - bphi0) ;
end