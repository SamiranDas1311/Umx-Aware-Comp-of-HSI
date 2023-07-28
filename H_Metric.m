function  Val=H_Metric(Im,alpha_1,alpha_2,N)
Im1= cat(3, Im, Im, Im);


[L,numLabels] = superpixels(Im1,N);


for i=1:numLabels
   e(i)=entropy(L(i));
   tv(i)=total_variation(L(i));
end
E=sum(e);
TV=sum(tv);
Val=alpha_1*E+alpha_2*TV;
end 




function TV = total_variation(y)
% y is the image
nbdims = 2;


% check number of channels in an image
if size(y,1)==1 || size(y,2)==1
    % we have one dimension
    nbdims = 1;
end


if size(y,1)>1 && size(y,2)>1 && size(y,3)>1
    % we have three dimensions
    nbdims = 3;
end


if nbdims==1
    TV = sum(abs(diff(y)));
    return;
end


% the total variation weight is 1
% weight_tv = ones(size(y));


g = gradient(y);
% compute using the l1 norm of the first order derivatives
TV = sum( abs(g),nbdims+1);


% TV = TV .* weight_tv;
TV = sum(TV(:));
end