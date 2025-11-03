function c = f_Consumption(aprime,a,z,alpha,delta,r)
% Action space: (aprime,a,z)
% State variables: (a,z)

w=(1-alpha)*((r+delta)/alpha)^(alpha/(alpha-1));

c=w*z+(1+r)*a-aprime; % Budget Constraint

end %end function