function F = Aiyagari1994_ReturnFn(aprime,a,z, alpha,delta,mu,r)
% Action space: (aprime,a,z)
% State variables: (a,z)

w=(1-alpha)*((r+delta)/alpha)^(alpha/(alpha-1));

c=w*z+(1+r)*a-aprime; % Budget Constraint

F=-Inf;
if c>0
    if mu==1
        F=log(c);
    else
        F=(c^(1-mu) -1)/(1-mu);
    end
end

end %end function