function HsRt = cvexTformToSRT(H)

R = H(1:2,1:2);
% Compute theta from mean of two possible arctangents
theta = mean([atan2(R(2),R(1)) atan2(-R(3),R(4))]);
% Compute scale from mean of two stable mean calculations
scale = mean(R([1 4])/cos(theta));
% Translation remains the same:
translation = H(3, 1:2);
% Reconstitute new s-R-t transform:
HsRt = [[scale*[cos(theta) -sin(theta); sin(theta) cos(theta)]; ...
  translation], [0 0 1]'];
% tformsRT = affine2d(HsRt);
% 
% imgBold = imwarp(imgB, tform, 'OutputView', imref2d(size(imgB)));
% imgBsRt = imwarp(imgB, tformsRT, 'OutputView', imref2d(size(imgB)));
% 
% figure(2), clf;
% imshowpair(imgBold,imgBsRt,'ColorChannels','red-cyan'), axis image;
% title('Color composite of affine and s-R-t transform outputs');
end