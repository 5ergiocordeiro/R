function plotparms()
	% Constantes
	mu = 4 * pi * 1e-7;		% H/m
	ep = 8.854 * 1e-12;		% F/m
	% Parâmetros fixos
	sigma = 0.1;			% H/m
	epr = 4.0;
	mur = 600.0;
	faixa = 0:110;			% log(1 Hz) a log(100 GHz)
	omega = 2 * pi * 10.^(faixa/10);
	k = omega .* (mur * mu * (epr * ep - j * sigma ./ omega)).^0.5;
	alpha = real(k);
	beta = imag(k);
	eta = j * mu * mur * omega ./ k;
	meta = abs(eta);
	teta = arg(eta);
	plot(faixa,alpha);
	plot(faixa,beta);
	plot(faixa,real(eta));		
	plot(faixa,imag(eta));
	plot(faixa,meta);			
	plot(faixa,teta);			
end
	