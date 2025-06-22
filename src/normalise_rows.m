function [A,b] = normalise_rows(A,b)
  normsA=sqrt(sum(A.^2,2));
  idx=normsA>0;
  A(idx,:)=bsxfun(@rdivide,A(idx,:),normsA(idx));
  b(idx)=b(idx)./normsA(idx);  
end

