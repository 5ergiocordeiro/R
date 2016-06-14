function Vmap = FEM_Tri_SolveP(V, ref, n, h)
	% Organiza os potenciais calculado e de referência em uma matriz
	Vmap = zeros(n+1,n+1);
	Vref = Vmap;
	global de
	for i = 1:n+1
		for j = 1:n+1
			ei = de(i, j);
			Vmap(i,j) = V(ei);
			Vref(i,j) = ref(i,j,h);
		end
	end
	% Plota a distribuição encontrada e a de referência
	ticks = linspace(1, n+1, 6);
	figure(1), contour(Vmap'), colorbar, title('Potencial (V)'); ax = gca; ax.XTick=ticks; ax.YTick=ticks; grid on;
	figure(2), contour(Vref'), colorbar, title('Potencial (V)'); ax = gca; ax.XTick=ticks; ax.YTick=ticks; grid on;
end