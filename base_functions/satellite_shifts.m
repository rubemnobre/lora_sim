function shifts = satellite_shifts(altitude, center_freq, init_elev, fs, n)
    % Valid for overhead passes and circular orbits, does not consider
    % earth rotation

    M = 5.972e24;   % Earth Mass
    RE = 6.378e6;   % Earth radius
    G = 6.6743e-11; % Gravitational constant
    C = 3e8;        % Speed of light

    rs = RE + altitude;         % Satellite orbit radius
    op = 2*pi*sqrt(rs^3/(G*M)); % Satellite orbital period
    ws = 2*pi/op;               % Satellite angular velocity

    initial_a = acos((RE - RE*sin(init_elev)^2 + sin(init_elev)*(rs^2 + RE^2*sin(init_elev)^2 - RE^2)^(1/2))/rs);
    final_a = initial_a - ws*n/fs;

    thetas = linspace(initial_a, final_a, n+1);
    ranges = sqrt(rs^2 + RE^2 - 2*rs*RE*cos(thetas));  % Slant range equation
    velocities = diff(ranges)*fs;

    beta = velocities/C;
    shifts = center_freq*sqrt((1 - beta)./(1 + beta)) - center_freq; % Relativistic doppler effect
end