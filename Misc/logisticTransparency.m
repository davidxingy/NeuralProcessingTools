function alphaData = logisticTransparency(pixelVals,cutoff,spread)
% alphaData = logisticTransparency(pixelVals,cutoff,spread)
% Function to make image/imagesc plots with transparent background for
% values close to 0 (with some anti-aliasing)

logFun = @(x,mu,sigma) 1./(1+exp(-1*(x-mu)/sigma));
alphaData = logFun(pixelVals,cutoff,spread);
alphaData(alphaData<0.1) = 0;
alphaData(alphaData>0.7) = 1;
