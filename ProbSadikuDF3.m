function Z0 = ProbSadikuDF3(h, nt)
% Resolve o exemplo 3.6 do livro "Numerical Techniques in Electromagnetic with MATLAB", de Matthew N. O. Sadiku, terceira edição, pp. 132 a 159}
	% Inicialização
	% ... constantes físicas:
	global epsilon0 c V0;
	[epsilon0, c] = deal(8.81e-12, 3e8);
	% ... parâmetros do problema:
	% ...... dimensões em cm
	[a, b, d, w, epsilonr] = deal(2.5, 2.5, 0.5, 1.0, 2.35);
	% ... posições fixas
	pos = num2cell(round([a b d w]/h));
	[nx, ny, jd, iw] = deal(pos{:});
	% Cálculos
	% ... calcula o potencial duas vezes, uma com e a outra sem o dielétrico
	V0 = ProbSadikuDF_CalcVI(nx, ny, jd, iw, 1, nt);
	Vd = ProbSadikuDF_CalcVI(nx, ny, jd, iw, epsilonr, nt);	
	% ... calcula o campo elétrico para cada distribuição de potencial
	E0 = Grad(V0, h * 0.01);
	Ed = Grad(Vd, h * 0.01);
	% ... plota os gráficos de V e E
	plotVE(V0, E0, 2);
	plotVE(Vd, Ed, 4);	
	% ... calcula a carga elétrica = capacitãncia, porque a diferença de potencial aplicada foi 1 Volt
	%out = round(([iw jd] + [nx ny]) / 2);
	out = [iw jd] + 3;
	out
	C0 = 4 * CalcQ(V0, out(1), out(2), jd, [1 1] * epsilon0);
	Cd = 4 * CalcQ(Vd, out(1), out(2), jd, [1 epsilonr] * epsilon0);
	% ... calcula a impedância da linha
	Z0 = 1 / (c * sqrt(C0 * Cd));
end

function V = ProbSadikuDF_CalcVI(nx, ny, jd, iw, epsilonr, nt)
% Calcula o potencial elétrico, por um processo iterativo, na região limitada por 0 <= i <= 'nx' e 0 <= j <= 'ny', com dielátrico de permissividade 'epsilonr' na posição (0 <= i <= 'nx', 0 <= j <= 'jd') e placa condutora na posição (0 <= i <= 'iw', j = 'jd').
	% Ajuste dos parâmetros
	% ... somar 1 porque os índices no MATLAB começam em 1, não em 0.
	% ... somar 1 porque os parâmetros passados são números de intervalos, não
	ajust = num2cell([nx, ny, jd, iw] + 2);
	[anx, any, ajd, aiw] = deal(ajust{:});
	% Inicializa o potencial
	global epsilon0;
	V = zeros(anx, any);
	V(1:aiw, ajd) = 1;
	epsilon = [1 epsilonr] * epsilon0;
	p = epsilon / (2 * sum(epsilon));
	% Recalcula o potencial 'nt' vezes
	for k=1:nt
		for i=1:(anx - 1)
			for j=1:(any - 1)
				if ( (j == ajd) && (i <= aiw) )
					% placa condutora; não é preciso recalcular V
				elseif (j == ajd)
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
	lap = zeros(anx-2,any-2);
	for i = 2:anx-1
		for j = 2:any-1
			lap(i,j) = V(i,j) - 0.25 * (V(i+1,j) + V(i-1,j) + V(i,j+1) + V(i,j-1));
		end
	end
	figure(7), imagesc(lap), colorbar;
end

function v = Grad(f, h)
% Calcula o gradiente em cada ponto do campo escalar 'f'. Usa diferenças progressivas de primeira ordem; nos limites superiores, repete o valor calculado para a penúltima linha/coluna.
	numxy = size(f);
	vx = zeros(numxy);
	vy = vx;
	for i = 1:(numxy(1) - 1)
		for j = 1:(numxy(2) - 1)
			vx(i,j) = (f(i+1,j) - f(i,j)) / h;
			vy(i,j) = (f(i,j+1) - f(i,j)) / h;
		end
		vx(i,end) = vx(i,end-1);
		vy(i,end) = vy(i,end-1);
	end
	vx(end,:) = vx(end-1,:);
	vy(end,:) = vy(end-1,:);
	v = {vx, vy};
end

function plotVE(V, E, fig)
% Plota o gráfico do potencial 'V' e do campo elétrico 'E'.
	% Equipotenciais
	figure(fig), contour(V'),  colorbar, title('Potencial elétrico (V)');
	% Vetores do campo elétrico
	% ... seleciona apenas 10 linhas e 10 colunas e amplifica 100 vezes o módulo, para melhor apresentação
	[Ex, Ey] = E{:};
	n = size(Ex);
	tam = min(n, [10 10]);
	idx = round(linspace(1, n(1), tam(1)));
	idy = round(linspace(1, n(2), tam(2)));
	selx = - 100 * Ex(idx,idy);
	sely = - 100 * Ey(idx,idy);
	figure(fig + 1), quiver(sely, selx), title('Campo elétrico');
end

function q = CalcQ(V, iout, jout, jd, epsilon)
% Calcula a carga por meio da lei de Gauss a partir do potencial 'V'
	% Ajuste dos parâmetros
	% ... somar 1 porque os índices no MATLAB começam em 1, não em 0.
	% ... somar 1 porque os parâmetros passados são números de intervalos, não índice dos pontos da grade
	ajust = num2cell([iout, jout, jd] + 2);
	[aiout, ajout, ajd] = deal(ajust{:});
	% Parâmetros para o cálculo
	% ... intervalos a considerar
	interval = [
		1, 1, ajout + 1, ajout + 1;
		1, 1, ajout, ajout;
		2, aiout, ajout + 1, ajout + 1;
		2, aiout, ajout, ajout;
		aiout + 1, aiout + 1, ajd + 1, ajout;
		aiout, aiout, ajd + 1, ajout;
		aiout + 1, aiout + 1, 2, ajd - 1;
		aiout, aiout, 2, ajd - 1;
		aiout + 1, aiout + 1, 1, ajd - 1;
		aiout, aiout, 1, ajd - 1;
		aiout + 1, aiout + 1, ajd, ajd;
		aiout, aiout, ajd, ajd;
		aiout + 1, aiout + 1, ajd, ajd;
		aiout, aiout, ajd, ajd;		
		];
	% ... pesos e permissividades em cada intervalo
	mult = epsilon([1 1 1 1 1 1 2 2 2 2 1 1 2 2]) .* [-0.5 0.5 -1 1 -1 1 -1 1 -0.5 0.5 -0.5 0.5 -0.5 0.5];
	% Cálculo pela lei de Gauss
	q = 0; 
	for k = 1:14
		q = q + sum(V(interval(k,1):interval(k,2), interval(k,3):interval(k,4))) * mult(k);
	end
end
