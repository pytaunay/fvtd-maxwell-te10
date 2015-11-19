%%% Function eh_te10
%%% returns the exact values of the TE10 mode for E and H
%%% based on geometry, time, and initial value for amplitude

function EH = eh_te10(x,z,t,a,A10,om,mu,eps,beta)

	kc2 = (pi/a)^2;

	E = zeros(3,1);
	E(2,1) = -1i*om*mu*pi/(kc2*a)*A10*sin(pi*x/a);
	
	H = zeros(3,1);
	H(1,1) = 1i*beta*pi/(kc2*a)*A10*sin(pi*x/a);
	H(3,1) = A10*cos(pi*x/a);

	EH = cat(1,E,H);
	EH = EH*exp(1i*om*t)*exp(-1i*beta*z);

end 


