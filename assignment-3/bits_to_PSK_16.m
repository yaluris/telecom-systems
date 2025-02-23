function X = bits_to_PSK_16(bit_seq)

dec_seq = bi2de(bit_seq,'left-msb');
X = zeros(length(dec_seq),2);

for i = 1:length(dec_seq)
    theta = (2*pi*dec_seq(i))/16;
    X(i,1) = cos(theta);
    X(i,2) = sin(theta);
end

end
