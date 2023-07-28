function [nEm,MixedNorm]=MixedNormEval(EigVal)
[NoBand]=size(EigVal);
for I=1:NoBand-1
   REA(i)=(EigVal(i)-EigVal(i+1))/(EigVal(i)+EigVal(i+1));
end
Th=0.01
[~,k_ub]=min(REA-Th);
EigValDiff=diff(EigVal);
MixedNorm=0;
for i=2:k_ub
    MixedNorm(i) = norm(EigValDiff(i:NoBand-1),1)/(norm(EigValDiff(i-1:NoBand-2,1),0.5)+0.01);
    i=i+1;
end
[~,nEm]=max(MixedNorm);
nEm=nEm+2;
end