function testcalcevk(letra)
%{
cd("C:/Users/ASUS.1/Dropbox/Particular/r");
cd "C:/Documents and Settings/Cliente/Meus documentos/Dropbox/Particular/r"
testcalcevk("a");
%}
%{
Calcula e plota o campo elétrico e o potencial elétrico para a distribuição de cargas indicada pela letra.
Gera arquivos de saída no formato LATEX e nome ("trab13figba" & letra & "1|2" & ".tex"). 
Uso:
	testcalcevk(letra)
	onde
		letra é a letra do exercício.
	
Testado com Octave 4.0 e ambiente gráfico "qt" (default).
Usar codificação ANSI no editor para que os caracteres acentuados apareçam corretamente.
%}
	warning("off");
	% Cores
	darkred = [139, 0, 0]/255;
	darkgreen = [0 139 0]/255;	
	violet = [111 0 255]/255;
	darkorange = [139 64 0]/255;
	darkyellow = [139 139 0]/255;
	somecolor = [64 139 0]/255;	
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
			end
	GEZ = GEX = GEY = GE = Gx = GV = [];
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
	figure(1); clf; plot2axis(Gx,GV,GE);
	figure(2); clf; plotvec(Gx,GEX,GEY,GEZ);
	
	
		[ax,pV,pE] = plotyy(Gx,GV,Gx,GE,'plot','plot');
		set(pV,'color',darkred,'linestyle','-'); set(pE,'color',darkgreen,'linestyle','-');
		set(ax(1),'ycolor',darkred); set(ax(2),'ycolor',darkgreen);
		set(ax(1),'ylim',limits([4 5])); set(ax(2),'ylim',limits([6 7]));
		xlabel(ax(1),'Posição x (m)'); ylabel(ax(1),'Potencial (V)'); ylabel(ax(2),'Módulo do campo elétrico (V/m))');
		legend("V","E");
	title("y = 0, z = 0");
%}
	% print -depslatex "-SX,Y" "trab13figbaa1.tex";
	figure(2); clf; hold on;
		plot(Gx,GEX,'color',darkyellow,'linestyle','--',
			Gx,GEY,'color',darkorange,'linestyle','--',
			Gx,GEZ,'color',somecolor,'linestyle','--');
		xlim([limits(1),limits(3)]);
		ylim([limits(6),limits(7)]);
		xlabel('Posição x (m)'); ylabel('Campo Elétrico (V/m)','color',darkgreen);
		legend("E_x","E_y","E_z");
		title("y = 0, z = 0");
	% print -depslatex "-SX,Y" "trab13figbaa2.tex";
%{
	figure(3); clf; 
	[ax,pEX,pEY] = plotyy(Gx,GEX,Gx,GEY,'plot','plot');
		set(pEX,'color',darkyellow,'linestyle','--'); set(pEY,'color',darkorange,'linestyle',':');
		set(ax(1),'ycolor',darkgreen);
		set(ax(1),'ylim',limits([6 7])); set(ax(2),'ylim',limits([6 7]));
		xlabel(ax(1),'Posição (m)'); ylabel(ax(1),'Campo Elétrico (V/m))');
		legend("E_x","E_y");
		nrows = size(Q)(1);
		for idx=1:nrows
			text(Q(idx,2), Q(idx,3),"*", 'horizontalalignment','center','verticalalignment','middle','color',violet);
			text(Q(idx,2), Q(idx,3),["Q_",num2str(idx)], 'horizontalalignment','center','verticalalignment','top','color','blue');
		end
%}