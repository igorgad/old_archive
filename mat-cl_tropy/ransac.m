function offset = ransac(pts, iterNum, thDist, thInlrRatio,  sampleNum)

ptNum = length(pts);
thInlr = round(thInlrRatio*ptNum);
inlrNum = zeros(1,iterNum);
ofst = zeros(1,iterNum);

for p = 1:iterNum
	% 1. fit using 2 random points
    sampleIdx = ceil(ptNum*rand(1,sampleNum));
	ptSample = pts(sampleIdx);
	d = mean(ptSample);
	
	% 2. count the inliers, if more than thInlr, refit; else iterate
	dist1 = pts-d;
	inlier1 = find(abs(dist1) < thDist);
	inlrNum(p) = length(inlier1);
	if length(inlier1) < thInlr, continue; end
	ofst(p) = mean(pts(:,inlier1));
    
end

% 3. choose the coef with the most inliers
[~,idx] = max(inlrNum);
offset = ofst(idx);

end