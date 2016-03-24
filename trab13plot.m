function trab13plot()
	special_str = [
		["Intensidade do campo eletrico em funcao de z"] ;
		["Intensidade do campo eletrico proximo a z = 0"];
		["Potencial eletrico em funcao de z"]
		];
	darkred = [139 0 0]/255;
	darkgreen = [0 139 0]/255;	
	% Problema 1
	% letra a
	R = 1;
	z = linspace(-10*R,10*R,100);
	r = R^2 + z.^2;
	E = z.* r.^(-3/2);
	m = 100/max(E);
	figure(1); clf; plot(z,E*m,'color',darkgreen);
		title(strtrim(special_str(1,:)));
		xlabel("z");
		ylabel("Emax");
	print -djpg -color trab13figbca.jpg;
	% letra b
	% geral	
	R = 1;
	z = linspace(-4*R,4*R,100);
	r = R^2 + z.^2;
	E = (sign(z) - (z.* r.^(-1/2)));
	m = 100/max(E);
	figure(1); clf; plot(z,E*m,'color',darkgreen);
		title(strtrim(special_str(1,:)));
		xlabel("z");
		ylabel("Emax");
		xlim([-4*R 4*R]);
	print -djpg -color trab13figbcb1.jpg;
	% detalhe
	R = 1;
	z = linspace(-0.5*R,0.5*R,100);
	r = R^2 + z.^2;
	E = (sign(z) - (z.* r.^(-1/2)));
	m = 100/max(E);
	figure(2); clf; plot(z,E*m,'color',darkgreen);
		title(strtrim(special_str(2,:)));
		xlabel("z");
		ylabel("Emax");
		xlim([-0.5*R 0.5*R]);
	print -djpg -color trab13figbcb2.jpg;
	% Problema 2
	L = 1;
	y = linspace(-4*L,4*L,100);
	u = y./L;
	f = sqrt(4 * u.^2 + 1);
	V = log((f.+1)./(f.-1));
	m = 100/max(V);
	figure(1); clf; plot(u,V*m,'color',darkred);
		title(strtrim(special_str(3,:)));
		xlabel("y/L");
		ylabel("Vmax");
	print -djpg -color trab13figbcd.jpg;