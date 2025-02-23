function [est_X,est_bit_seq] = detect_PSK_16(Y)

est_X = zeros(length(Y),2);
est_bit_seq = zeros(length(Y),4);

angle = zeros(1,16);
for m = 0:15
    angle(m+1) = (2*pi*m)/16;
end

for i = 1:length(Y)
    theta = atan2(Y(i,2),Y(i,1));
    angularDifferences = abs(mod(angle-theta+pi,2*pi)-pi);
    [~,index] = min(angularDifferences);
    est_X(i,1) = cos(angle(index));
    est_X(i,2) = sin(angle(index));
    est_bit_seq(i,:) = de2bi(index-1,4,'left-msb');
end

end
