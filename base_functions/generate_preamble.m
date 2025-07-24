function s = generate_preamble(SF, B, OSR, NUP, SW)
    Ns = 2^SF*OSR;
    N = round(Ns*(NUP+4.25));
    s = zeros(1, N);
    s0 = generate_lora_symbol(0, SF, B, OSR, false, 0);
    s1 = generate_lora_symbol(double(bitshift(SW, -4)*8), SF, B, OSR, false, 0);
    s2 = generate_lora_symbol(double(bitand(SW, 0xF)*8), SF, B, OSR, false, 0);
    for i=1:NUP
        s(1 + (i-1)*Ns: i*Ns) = s0;
    end
    s(1 + NUP*Ns:(NUP+1)*Ns) = s1;
    s(1 + (NUP+1)*Ns:(NUP+2)*Ns) = s2;
    s(1 + (NUP+2)*Ns:(NUP+3)*Ns) = conj(s0);
    s(1 + (NUP+3)*Ns:(NUP+4)*Ns) = conj(s0);
    s(1 + (NUP+4)*Ns:(NUP+4)*Ns + Ns/4) = conj(s0(1:Ns/4));
end