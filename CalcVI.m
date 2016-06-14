function V = CalcVI(nx, ny, jd, iw, epsilonr, nt)
% Calcula o potencial elétrico, por um processo iterativo, na região limitada por 0 <= i <= 'nx' e 0 <= j <= 'ny', com dielátrico de permissividade 'epsilonr' na posição (0 <= i <= 'nx', 0 <= j <= 'jd') e placa condutora na posição (0 <= i <= 'iw', j = 'jd').
	% Inicializa o potencial
	V = zeros(nx, ny);
	V(1:iw, jd) = 1;
	epsilon = [1 epsilonr];
	p = 2 * epsilon / (sum(epsilon));
	% Recalcula o potencial 'nt' vezes
	for k=1:nt
		for i=1:(nx - 1)
			for j=1:(ny - 1)
				if ( (j == jd) && (i <= iw) )
					% placa condutora; não é preciso recalcular V
				elseif (j == jd)
					% interface entre o substrato e o ar
					V(i,j) = 0.25 * (V(i+1,j) + V(i-1,j) + p(1)* V(i,j+1) + p(2)* V(i,j-1));
				elseif ((i == 1) && (j == 1))
					% origem
					V(i,j) = 0.25 * (2 * V(i+1,j) + 2 * V(i,j+1));
				elseif (i == 1)
					% eixo Y
					V(i,j) = 0.25 * (2 * V(i+1,j) + V(i,j+1) + V(i,j-1));
				elseif (j == 1)
					% eixo X
					V(i,j) = 0.25 * (V(i+1,j) + V(i-1,j) + 2 * V(i,j+1));
				else
					% ar ou substrato
					V(i,j) = 0.25 * (V(i+1,j) + V(i-1,j) + V(i,j+1) + V(i,j-1));
				end
			end
		end
		% remover os comentários da linha abaixo para plotagem dos resultados parciais
		%figure(1), imagesc(flipud(V')), colorbar, title([num2str(k), '/', num2str(nt)]);	drawnow;
	end
end

