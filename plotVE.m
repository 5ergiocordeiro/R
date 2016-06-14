function plotVE(V, E, fig)
% Plota o gráfico do potencial 'V' e do campo elétrico 'E' na figura 'fig'.
	feature('DefaultCharacterSet', 'UTF8');
	% Equipotenciais
	n = size(V, 1);
	idx = linspace(1, n, n);
	[X, Y] = meshgrid(idx, idx);
	hf = figure(fig);
	contour(X, Y, V'), colorbar;
	% Vetores do campo elétrico
	% ... seleciona apenas 20 linhas e 20 colunas, para melhor apresentação
	[Ex, Ey] = E{:};
	n = size(Ex);
	tam = min(n, [20 20]);
	idx = round(linspace(1, n(1), tam(1)));
	idy = round(linspace(1, n(2), tam(2)));
	[X, Y] = meshgrid(idx, idy);
	selx = - Ex(idx, idy);
	sely = - Ey(idx, idy);
	figure(fig), hold on, quiver(Y, X, selx, sely), title('Campo elétrico');
	hold off;
	saveas(hf, strcat('plotVE', num2str(fig), '.jpg'));
end
