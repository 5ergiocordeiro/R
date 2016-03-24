function [V,calc]=initv(Rext,Rint,Lext,Lint,rlim,zlim,Vext,Vint,RhoNint,RhoNext,ZNint,ZNext,Nrho,Nz);
% Inicializa o potencial 'V' para o cálculo iterativo.
% Retorna a grade inicializada e a posição dos cilindros nela.
% Chamada pela função 'calcv'.
	% Cria a grade
	V = zeros(Nrho,Nz);
	% Calcula a aproximação inicial.
	% Teoricamente, a aproximação inicial pode ser qualquer; no etantom para tentar auxiliar a convergência, foram escolhidas as condições seguintes:
	% ... Vint dentro do cilindro interno
	V(1:RhoNint,ZNint(1):ZNint(2)) = Vint;
	% ... queda exponencial até Vext no cilindro externo
	a = (Vint - Vext) / (RhoNint - RhoNext);
	for i = RhoNint:RhoNext
		V(i,ZNext(1):ZNext(2)) = Vint * log(RhoNext/i);
	end
	% Onde não calcular o Laplaciano
	calc = ones(Nrho,Nz);
	calc(RhoNext,ZNext(1):ZNext(2)) = 0 ;		% no cilindro externo
	calc(RhoNint,ZNint(1):ZNint(2)) = 0;		% no cilindro interno
end
