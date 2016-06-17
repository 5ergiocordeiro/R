function Z0 = ProbSadikuDF(h, nt)
	e0 = 8.81e-12; a = 2.5; b = 2.5; d = 0.5; w = 1.0; er = 2.35;
	nx = round(a/h); ny = round(b/h); jd = round(d/h); iw = round(w/h);
	V0 = ProbSadikuDF_CalcVI(nx, ny, jd, iw, 1, nt);
	E0 = Grad(V0, h);
	plotVE(V0, E0,2);
	V1 = ProbSadikuDF_CalcVI(nx, ny, jd, iw, er, nt);	
	E1 = Grad(V1, h);
	plotVE(V1, E1,4);	
	iout = round((iw + nx)/2); jout = round((jd + ny)/2);
	%iout = iw + 3; jout = jd + 3;
	C0 = 4 * CalcQ(V0, iout, jout, jd, [e0, e0]);
	C1 = 4 * CalcQ(V1, iout, jout, jd, [e0, er*e0]);
	Z0 = 1 / (3e8 * sqrt(C0 * C1));
end

function V = ProbSadikuDF_CalcVI(nx, ny, jd, iw, er, nt)
	V = zeros(nx+2,ny+2);
	V(2:iw+2,jd+2) = 1;
	p1 = 2 / (1 + er);
	p2 = 2 * er / (1 + er);
	for k=1:nt
		for i=0:nx - 1
			for j=0:ny - 1
				if ( (j == jd) && (i <= iw) )
				elseif (j == jd)
					V(i+2,j+2) = 0.25 * (V(i+3,j+2) + V(i+1,j+2) + p1*V(i+2,j+3) + p2*V(i+2,j+1));
				elseif ((i == 0) && (j == 0))
					% origem
					V(i+2,j+2) = 0.25 * (2 * V(i+3,j+2) + 2 * V(i+2,j+3));
				elseif (i == 0)
					% eixo Y
					V(i+2,j+2) = 0.25 * (2 * V(i+3,j+2) + V(i+2,j+3) + V(i+2,j+1));
				elseif (j == 0)
					% eixo X
					V(i+2,j+2) = 0.25 * (V(i+3,j+2) + V(i+1,j+2) + 2 * V(i+2,j+3));
				else
					% ar ou substrato
					V(i+2,j+2) = 0.25 * (V(i+3,j+2) + V(i+1,j+2) + V(i+2,j+3) + V(i+2,j+1));
				end
			end
		end
		% remover os comentários da linha abaixo para plotagem dos resultados parciais
		%figure(1), imagesc(flipud(V')), colorbar, title([num2str(k), '/', num2str(nt)]);	drawnow;
	end
end

function v = Grad(f, h)
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
	figure(fig), contour(V(2:end,2:end)'),  colorbar, title('Potencial elétrico (V)');
	[Ex, Ey] = E{:};
	n = size(Ex);
	tam = min(n, [10 10]);
	idx = round(linspace(1, n(1), tam(1)));
	idy = round(linspace(1, n(2), tam(2)));
	selx = - Ex(idx,idy);
	sely = - Ey(idx,idy);
	figure(fig + 1), quiver(sely, selx), title('Campo elétrico');
end


function q = CalcQ(V, iout, jout, jd, e)
	for k = 1:2
		asum = e(1) * sum(V(3:iout+1, jout+2)) + e(1) * V(2,jout+2) / 2 + e(2) * V(iout+2, 2) / 2;
		for j = 2:jout+2
			if (j < jd)
				asum = asum + e(2) * V(iout + 2, j + 2);
			else
				asum = asum + e(1) * V(iout + 2, j + 2);
			end
		end
		if (k == 1)
			sv = asum;
		end
		iout = iout - 1;
		jout = jout - 1;
	end
	asum = asum + 2 * e(1) * V(iout + 2, jout + 2);
	q = sv - asum;
end

%{
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
%}