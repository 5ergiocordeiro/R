function plotvec(X,VX,VY,VZ,limits,texts,colors)
%{
%}
	plot(X,VX,'color',colors(1,:),'linestyle','--',
			X,VY,'color',colors(2,:),'linestyle','--',
			X,VZ,'color',colors(3,:),'linestyle','--');
		xlim([limits(1),limits(2)]);
		ylim([limits(3),limits(4)]);
		xlabel(strtrim(texts(1,:)));
		ylabel(strtrim(texts(2,:)),'color',colors(4,:));
		legend(strtrim(texts(3,:)),strtrim(texts(4,:)),strtrim(texts(5,:)));
		title(strtrim(texts(6,:)));
