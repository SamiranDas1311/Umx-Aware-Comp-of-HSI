function [X_Rec]=Recons(Fac, Lib,LibIdx)
A_Est=Fac{0}.*Fac{1].*Fac{2};
[nRow,nCol,nEm=size(A_Est);
Em_Lib=Lib(LibIdx,:);
A_2d=reshape(A_Est,[nRow*nCol,nEm]);
X_Rec=A_2d*Em_Lib;
end