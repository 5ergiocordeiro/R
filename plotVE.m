function plotVE(V, E, fig)
% Plota o gr�fico do potencial 'V' e do campo el�trico 'E' na figura 'fig'.
	feature('DefaultCharacterSet', 'UTF8');
	% Equipotenciais
	n = size(V, 1);
	idx = linspace(1, n, n);
	[X, Y] = meshgrid(idx, idx);
	hf = figure(fig);
	contour(X, Y, V'), colorbar;
	% Vetores do campo el�trico
	% ... seleciona apenas 20 linhas e 20 colunas, para melhor apresenta��o
	[Ex, Ey] = E{:};
	n = size(Ex);
	tam = min(n, [20 20]);
	idx = round(linspace(1, n(1), tam(1)));
	idy = round(linspace(1, n(2), tam(2)));
	[X, Y] = meshgrid(idx, idy);
	selx = - Ex(idx, idy);
	sely = - Ey(idx, idy);
	figure(fig), hold on, quiver(Y, X, selx, sely), title('Campo el�trico');
	hold off;
	saveas(hf, strcat('plotVE', num2str(fig), '.jpg'));
end
