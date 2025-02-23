function X = bits_to_4PAM(b)

X = zeros(1,length(b));
for i = 1:length(b)
    if (b(i,1) == 0) && (b(i,2) == 0)
        X(i) = 3;
    elseif (b(i,1) == 0) && (b(i,2) == 1)
        X(i) = 1;
    elseif (b(i,1) == 1) && (b(i,2) == 1)
        X(i) = -1;
    else
        X(i) = -3;
    end
end
