%script for testing cosine distance


foo = 0:.01:1;
tau = .5:.1:1;


figure; 
for i = 1:length(tau)
    W = exp(-(1./(foo+.01)-.9901)/tau(i));
    plot(foo,W,'Color',[0,0,tau(i)/max(tau)]);
    hold on;
end
hold off;


figure; 
for i = 1:length(tau)
    W = exp(-abs((1./(foo+.01)-.9901).^(1.1))/tau(i));
    plot(foo,W,'Color',[0,0,tau(i)/max(tau)]);
    hold on;
end
hold off;


figure; 
for i = 1:length(tau)
    W = exp(-abs((1./(foo+.01)-.9901).^(1.3))/tau(i));
    plot(foo,W,'Color',[0,0,tau(i)/max(tau)]);
    hold on;
end
hold off;