function plote(Erho,Ez,rlim,zlim,Nrho,Nz)
% Plota e grava os gráficos que ilustram o campo elétrico 'E' para a grade dada ('rlim','zlim','Nrho','Nz').
% Chamada pela função principal ('probcil').
	Eabs = (Erho.^2 + Ez.^2).^0.5;
	txV = linspace (0,rlim,Nrho);
	tyV = linspace (0,zlim,Nz);
	txE = txV(1:50:Nrho); 
	tyE = tyV(1:50:Nz);
	Ey = Erho(1:50:Nrho,1:50:Nz);
	Ex = Ez(1:50:Nrho,1:50:Nz);
	% Contorno
	figure(1); clf; hold on;
		contour(txV,tyV,Eabs);
		quiver(txE,tyE,Ex,Ey);
		ab = colorbar();
		set(ab,'title','V/m');
		print -djpg probcile1.jpg;
	% Heatmap
	figure(2); clf; hold on;
		pcolor(txV,tyV,Eabs);
		shading interp;
		ab = colorbar();
		set(ab,'title','V/m');
		print -djpg probcile2.jpg;
	% 3D
	figure(3); clf; hold on;
		meshc(tyV,txV,Eabs);
		ab = colorbar();
		set(ab,'title','V/m');		
		print -djpg probcile3.jpg;
	% 2D
	% Intensidade do campo para 4 valores diferentes de z
	figure(4); clf; hold on;
		plot(txV,Eabs(:,1),'r');
		plot(txV,Eabs(:,floor(ZNint/2)),'m');
		plot(txV,Eabs(:,ZNint),'k');
		plot(txV,Eabs(:,ZNext),'g');
		legend('meio dos cilindros','3/4 do cilindro interno','ponta do cilindro interno','ponta do cilindro externo');
		xlabel('r'); ylabel('|E|');
		print -djpg probcile4.jpg;
	% Intensidade do campo para 2 valores diferentes de rho
	figure(5); clf; hold on;
		plot(tyV,Eabs(floor(RhoNint/2),:),'c');
		plot(tyV,Eabs(floor((RhoNint+RhoNext)/2),:),'b');
		legend('dentro do cilindro interno','entre os cilindros');
		xlabel('z'); ylabel('|E|');
		print -djpg probcile5.jpg;
end
