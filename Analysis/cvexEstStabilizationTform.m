function H = cvexEstStabilizationTform(imgA, imgB, pointsA, pointsB)

imgA = rgb2gray(imgA);
imgB = rgb2gray(imgB);

% Extract FREAK descriptors for the corners
[featuresA, pointsA] = extractFeatures(imgA, pointsA);
[featuresB, pointsB] = extractFeatures(imgB, pointsB);

indexPairs = matchFeatures(featuresA, featuresB);
pointsA = pointsA(indexPairs(:, 1), :);
pointsB = pointsB(indexPairs(:, 2), :);

[tform, pointsBm, pointsAm] = estimateGeometricTransform(pointsB, pointsA, 'affine');

H = tform.T;
end