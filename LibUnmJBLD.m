function [XCov,YCov,DJbld,LibIdxJbld]=LibUnmJBLD(X,Lib,n)
[NoPix,NoBand]=size(X);
[NoLibEle,NoBand]=size(Lib);
XMean=mean(X);
XNorm=X-repmat(XMean,NoPix,1);
XCov=cov(X);

for i=1:NoLibEle
Y=[X;Lib(i,:)];
YMean=mean(Y);
YNorm=Y-repmat(YMean,NoPix+1,1);
YCov=cov(Y);
DJbld(i,1)=CovarianceSimMetric(XCov,YCov,5);
i=i+1;
end
%% Find the indices having the lowest covariance similarity distances
[r,idx]=sort(DJbld,'descend');
LibIdxJbld=idx(1:n,1);
end 