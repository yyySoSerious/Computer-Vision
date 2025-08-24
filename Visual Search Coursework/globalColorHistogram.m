function desc=globalColorHistogram(img, ndiv)
    qimg = double(img)./256; %r, g, b in [0, 1)
    qimg = floor(qimg.*ndiv); %r, g, b in [0, ndiv-1]

    bin = qimg(:, :, 1) * ndiv.^2 + qimg(:, :, 2) * ndiv + qimg(:, :, 3); %bin index in [0, ndiv^3 -1] (base ndiv)
    bin = reshape(bin, 1, []);
    desc = histogram(bin, ndiv.^3).Values;
    desc = desc ./ sum(desc);
end