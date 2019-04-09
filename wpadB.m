function [a,f1,b1]=wpad(a);
a=a';
nn=floor(log(length(a))/log(2));

mm=length(a);
axl=(2^(nn+1))-mm;
axm=floor(axl/2);

if (mod(axl,2)==0) 
   f1(1:axm)=a(1);
   b1(1:axm)=a(mm);
   a=[f1 a b1];
else
   f1(1:axm)=a(1);
   b1(1:axm+1)=a(mm);
   a=[f1 a b1];
end