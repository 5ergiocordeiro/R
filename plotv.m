function plotv(V,RhoNext,RhoNint,ZNext,ZNint,rlim,zlim,Nrho,Nz)
% Plota e grava os gráficos que ilustram o potencial 'V' para a geometria dada ('RhoNext','RhoNint','ZNext','ZNint','rlim','zlim','Nrho','Nz').
% Chamada pela função principal ('probcil').
	txV = linspace (0,rlim,Nrho);
	tyV = linspace (0,zlim,Nz);
	% Contorno
	figure(1); clf; hold on;
		contour(txV,tyV,V);
		ab = colorbar();
		set(ab,'title','V');
		print -djpg probcilv1.jpg;
	% Heatmap
	figure(2); clf; hold on;
		pcolor(txV,tyV,V);
		shading interp;
		ab = colorbar();
		set(ab,'title','V');
		print -djpg probcilv2.jpg;
	% 3D
	figure(3); clf; hold on;
		meshc(tyV,txV,V);
		ab = colorbar();
		set(ab,'title','V');
		print -djpg probcilv3.jpg;
	% 2D
	% Potencial para 4 valores diferentes de z
	figure(4); clf; hold on;
		plot(txV,V(:,floor(Nz/2)),'r');
		plot(txV,V(:,floor(Nz/2 + ZNint(2)/4)),'m');
		plot(txV,V(:,ZNint(2)),'k');
		plot(txV,V(:,ZNext(2)),'g');
		legend('meio dos cilindros','3/4 do cilindro interno','ponta do cilindro interno','ponta do cilindro externo');
		xlabel('r'); ylabel('V');
		print -djpg probcilv4.jpg;
	% Potencial para 2 valores diferentes de rho
	figure(5); clf; hold on;
		plot(tyV,V(floor(RhoNint/2),:),'c');
		plot(tyV,V(floor((RhoNint+RhoNext)/2),:),'b');
		legend('dentro do cilindro interno','entre os cilindros');
		xlabel('z'); ylabel('V');
		print -djpg probcilv5.jpg;
end
