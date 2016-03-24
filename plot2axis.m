function plot2axis(X,Y1,Y2,limits,texts,colors)
%{
%}
	[ax,Y1,Y2] = plotyy(X,Y1,X,Y2,'plot','plot');
		set(Y1,'color',colors(1,:),'linestyle','-');
		set(Y2,'color',colors(2,:),'linestyle','-');
		set(ax(1),'ycolor',colors(1,:));
		set(ax(2),'ycolor',colors(2,:));
		set(ax(1),'xlim',limits([1 2]));
		set(ax(1),'ylim',limits([3 4]));
		set(ax(2),'ylim',limits([5 6]));
		xlabel(ax(1),strtrim(texts(1,:)));
		ylabel(ax(1),strtrim(texts(2,:)));
		ylabel(ax(2),strtrim(texts(3,:)));
		legend(strtrim(texts(4,:)),strtrim(texts(5,:)));
		title(strtrim(texts(6,:)));
