function v = Grad(f, h)
% Calcula o gradiente do campo escalar 'f' em cada ponto da grade de largura 'h'. Usa diferenças progressivas de primeira ordem; nos limites superiores, repete o valor calculado para a penúltima linha/coluna.
	numxy = size(f);
	vx = zeros(numxy);
	vy = vx;
	vx(1:(numxy(1) - 1), 1:(numxy(2) - 1)) = (f(2:(numxy(1)), 1:(numxy(2) - 1)) - f(1:(numxy(1) - 1), 1:(numxy(2) - 1))) / h;
	vy(1:(numxy(1) - 1), 1:(numxy(2) - 1)) = (f(1:(numxy(1) - 1), 2:(numxy(2))) - f(1:(numxy(1) - 1), 1:(numxy(2) - 1))) / h;
	vx(1:(numxy(1) - 1), end) = vx(1:(numxy(1) - 1), end - 1);
	vy(1:(numxy(1) - 1), end) = vy(1:(numxy(1) - 1), end - 1);
	vx(end,:) = vx(end-1,:);
	vy(end,:) = vy(end-1,:);
	v = {vx, vy};
end
