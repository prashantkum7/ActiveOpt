function [out,f] = svmdecision (Xnew,svmstruct)


sv = svmstruct.SupportVectors;
alphaHat=svmstruct.Alpha;
bias = svmstruct.Bias;
kfun = svmstruct.KernelFunction;
kfunargs = svmstruct.KernelFunctionArgs;

f = (feval(kfun,sv,Xnew,kfunargs{:})'*alphaHat(:)) + bias;
out=sign(f);
out(out==0) = 1;
end