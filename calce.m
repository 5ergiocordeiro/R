function [Erho,Ez]=calce(V,rlim,zlim)
% Calcula o campo elétrico a partir do potencial, para a grade indicada ('rlim','zlim').
% Considera 'rlim' e 'zlim' em mm.
% Chamada pela função principal ('probcil.m').
	[nrow,ncol] = size(V);
	Erho = Ez = zeros(nrow,ncol);
	mrho = - 500 * nrow / rlim;
	mz = - 500 * ncol / zlim;
	Erho(2:nrow-1,2:ncol-1) = mrho * (V(3:nrow,2:ncol-1) - V(1:nrow-2,2:ncol-1));
	Ez(2:nrow-1,2:ncol-1) = mz * (V(2:nrow-1,3:ncol) - V(2:nrow-1,1:ncol-2));
end