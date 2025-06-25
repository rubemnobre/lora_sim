function a = e_to_a(e, H)
    syms alp elev;
    RE = 6.378e6;   % Earth radius
    d = sqrt(RE^2 + (RE+H)^2 - 2*RE*(RE+H)*cos(alp));
    sol = vpasolve([sin(alp)/d == sin(elev + pi/2)/(RE+H), elev==e]);
    a = double(sol.alp);
end