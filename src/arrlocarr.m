function [ind] = arrlocarr(A,B)
%arrlocarr finds index of common elements array A and B
%   Taken from https://www.mathworks.com/matlabcentral/answers/87731-find-elements-of-an-array-in-another-array
    [ii,jj]=unique(B,'stable');
    n=numel(A);
    ind = zeros(n,1);
    for k = 1:n
         ind(k)=find(ii == A(k));
    end
end

