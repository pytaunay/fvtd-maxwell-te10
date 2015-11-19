%%% Returns beta for the TE10 mode beta  = sqrt(k^2 - kc^2) = sqrt( om^2*mu*eps - (pi/a)^2)
function res = beta_te10(om,mu,eps,a)

	kc = pi/a;
	k = om*sqrt(mu*eps);
	
	res = sqrt(k^2-kc^2);

end


