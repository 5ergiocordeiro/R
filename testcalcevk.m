function testcalcevk(letra)
%{
cd("C:/Users/ASUS.1/Dropbox/Particular/r");
cd "C:/Documents and Settings/Cliente/Meus documentos/Dropbox/Particular/r"
testcalcevk("a");
%}
%{
Testa a função calcevk.
Calcula e plota o campo elétrico e o potencial elétrico para a distribuição de cargas indicada.
Gera arquivos de saída no formato LATEX e nome "trab14figba" letra "[1-8]" ".tex".
	Figuras:
	1: V e |E| ao longo do eixo x
	2: Componentes de E ao longo do eixo x
	3: V e |E| ao longo do eixo y
	4: Componentes de E ao longo do eixo y
	5: V no plano xy (3D)
	6: V no plano xy (heatmap)
	7: V e E no plano xy (contorno e vetores)
	8: E no espaço (vetores em 3D)
Uso:
	testcalcevk(letra)
	onde
		letra indica a distribuição de cargas a ser usada (ver comentário 'Distribuições de carga').
Utiliza funções de apoio: plot2axis, plotvec e plotheatmap.
Usar codificação ANSI no editor para que os caracteres acentuados apareçam corretamente.
Testado com Octave 4.0 e ambiente gráfico "qt" (default).
%}
	warning("off");				% devido ao excesso de zelo do interpretador
	% Definições
	% nome do arquivo de saída
	fname = "trab14figbaa";
	% resoluções
	gridsize = [50 20 3];
	% cores
	darkred = [139 0 0]/255;
	darkgreen = [0 139 0]/255;	
	violet = [111 0 255]/255;
	darkorange = [139 64 0]/255;
	darkyellow = [139 139 0]/255;
	somecolor = [64 139 0]/255;
	% textos fixos
	textsfig1 = [
		"Posição x (m)";"Potencial elétrico (V)"; "Módulo do campo elétrico (V/m))";
		"V"; "E"; "Eixo x (y = 0; z = 0)"
		];
	textsfig2 = [
		"Posição x (m)"; "Campo elétrico (V/m)";
		"E_x"; "E_y"; "E_z"; "Eixo x (y = 0; z = 0)"
		];		
	textsfig3 = [
		"Posição y (m)";"Potencial elétrico (V)"; "Módulo do campo elétrico (V/m))";
		"V"; "E"; "Eixo y (x = 0; z = 0)"
		];
	textsfig4 = [
		"Posição y (m)"; "Campo elétrico (V/m)";
		"E_x"; "E_y"; "E_z"; "Eixo y (x = 0; z = 0)"
		];		
	textsfig6 = [
		"Posição x (m)"; "Posição y (m)"; "Potencial no plano xy (z = 0)"; "V"
		];		
	textsfig7 = [
		"Posição x (m)"; "Posição y (m)"; "Potencial e campo elétrico no plano xy (z = 0)"; "V"; "E"
		];		
	textsfig8 = [
		"Posição x (m)"; "Posição y (m)"; "Campo elétrico";
		];		
	colorsfig1 = [darkred;darkgreen];
	colorsfig2 = [darkyellow;darkorange;somecolor;darkgreen];
	% distribuições de carga (indicadas por letra) e respectivos limites
	switch letra
		case "a"
			Q = [140e-3 -60e-3 0 0 ; -140e-3 60e-3 0 0 ];
			limits = [-0.3 0.001 0.3 -100 100 -2500 2500];
		case "b"
			Q = [140e-3 -60e-3 0 0 ; 140e-3 60e-3 0 0 ];
			limits = [-0.3 0.001 0.3 -100 100 -2500 2500];
		case "c"
			Q = [100e-3 -0.03 -0.03 0; -100e-3 -0.03 0.03 0; -100e-3 0.03 -0.03 0; 100e-3 0.03 0.03 0];
			limits = [-0.2 0.002 0.2 -0.1 0.1 -2500 2500];
		case "d"
			Q = [100e-3 -0.03 -0.03 0];
			limits = [-0.2 0.002 0.2 -0.1 0.1 -2500 2500];
			end
	% Fim das definições
%{
	% Figuras 1 e 2
	% Calcula os valores de V e E no intervalo dado por limits
	GEZ = GEX = GEY = GE = Gy = Gx = GV = [];
	for x = limits(1):limits(2):limits(3)
		[erro,V,E] = calcevk( [x 0 0], Q ) ;
		if erro == 0
			Gx = [Gx,x];
			GV = [GV,V];
			GEX = [GEX,E(1)];
			GEY = [GEY,E(2)];
			GEZ = [GEZ,E(3)];
			GE = [GE,norm((E),2)];
		end
	end
	figure(1); clf; plot2axis(Gx,GV,GE,limits([1 3:7]),textsfig1,colorsfig1);
	% print -depslatex "-SX,Y" fname "1.tex";
	figure(2); clf; plotvec(Gx,GEX,GEY,GEZ,limits([1 3 6 7]),textsfig2,colorsfig2);
	% print -depslatex "-SX,Y" fname "2.tex";
	% Figuras 3 e 4
	% Calcula os valores de V e E no intervalo dado por limits
	GEZ = GEX = GEY = GE = Gy = Gx = GV = [];
	for y = limits(1):limits(2):limits(3)
		[erro,V,E] = calcevk( [0 y 0], Q ) ;
		if erro == 0
			Gy = [Gy,y];
			GV = [GV,V];
			GEX = [GEX,E(1)];
			GEY = [GEY,E(2)];
			GEZ = [GEZ,E(3)];
			GE = [GE,norm((E),2)];
		end
	end
	figure(3); clf; plot2axis(Gy,GV,GE,limits([1 3:7]),textsfig3,colorsfig1);
	% print -depslatex "-SX,Y" fname "3.tex";
	figure(4); clf; plotvec(Gy,GEX,GEY,GEZ,limits([1 3 6 7]),textsfig4,colorsfig2);
	% print -depslatex "-SX,Y" fname "4.tex";
%}
	% Figuras 5 a 8
	% Calcula os valores de V no intervalo dado por limits, mas com menor resolução que no caso das figuras anteriores
	GV = zeros(50,50);
	txV = tyV = linspace (limits(1), limits(3), 50)';
	for idx = 1:50
		for idy = 1:50
			[erro,V,E] = calcevk( [txV(idx) tyV(idy) 0], Q ) ;
			GV(idx,idy) = V;
		end
	end
	% Calcula os valores de E no intervalo dado por limits, mas com resolução ainda menor
	GE = GEX = GEY = GEZ = zeros(20,20);
	txE = tyE = linspace (limits(1), limits(3), 20)';
	for idx = 1:20
		for idy = 1:20
			[erro,V,E] = calcevk( [txE(idx) tyE(idy) 0], Q ) ;
			GEX(idx,idy) = E(1);
			GEY(idx,idy) = E(2);
			GEZ(idx,idy) = E(3);
			GE(idx,idy) = norm((E),2);
		end
	end
	% figure(5); clf; mesh(txV, tyV, GV);
	% print -depslatex "-SX,Y" fname "5.tex";
	% figure(6); clf; plotheatmap(txV,tyV,GV.',textsfig6,'hot','blur','');
	% print -depslatex "-SX,Y" fname "6.tex";
	figure(7); clf; hold on;
		plotheatmap(txV,tyV,GV.',textsfig7,'hot','contour','');
		quiver (txE,tyE,GEX.',GEY.','color',darkgreen);
		title(strtrim(textsfig7(3,:)));
		xlim([-0.2 0.2]); ylim([-0.2 0.2]);
	% print -depslatex "-SX,Y" fname "7.tex";
%{
	% Calcula os valores de V e E no intervalo dado por limits, mas com resolução ainda menor
	GV = GE = GEX = GEY = GEZ = zeros(20,20,3);
	tx = ty = linspace (limits(1), limits(3), 20)';
	tz = linspace (limits(1), limits(3), 3)';
	for idz = 1:3
		for idx = 1:20
			for idy = 1:20
				[erro,V,E] = calcevk( [tx(idx) ty(idy) tz(idz)], Q ) ;
				GV(idx,idy,idz) = V;
				GEX(idx,idy,idz) = E(1);
				GEY(idx,idy,idz) = E(2);
				GEZ(idx,idy,idz) = E(3);
				GE(idx,idy,idz) = norm((E),2);
			end
		end
	end
	figure(10); clf;
		quiver3 (tx, ty, tz , GEX,GEY,GEZ,'color',darkgreen);
		title(strtrim(textsfig9(3,:)));
		xlim([-0.2 0.2]); ylim([-0.2 0.2]); 
%}
