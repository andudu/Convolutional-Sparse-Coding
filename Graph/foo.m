% try and plot the potential l1 weights

l1 = L{1}.ind1(1):L{1}.ind2(1);
r1 = L{1}.ind1(2):L{1}.ind2(2);
x = zeros(length(l1),length(r1),16);
y = zeros(L{end}.ind2(1),L{end}.ind2(2));

for i = 1:16
    for j = 1:4
        if i == 2
            if j <=3
                temp = abs(L{i}.phi(:,j));
                l1 = L{i}.ind1(1):L{i}.ind2(1);
                r1 = L{i}.ind1(2):L{i}.ind2(2);
                x(:,:,i) = reshape(temp,length(l1),length(r1))+x(:,:,i);
                y(l1,r1) = y(l1,r1) + reshape(temp,length(l1),length(r1));
            end
        else
            temp = abs(L{i}.phi(:,j));
            l1 = L{i}.ind1(1):L{i}.ind2(1);
            r1 = L{i}.ind1(2):L{i}.ind2(2);
            x(:,:,i) = reshape(temp,length(l1),length(r1))+x(:,:,i);
            y(l1,r1) = y(l1,r1) + reshape(temp,length(l1),length(r1));
        end
    end
end
square_plot(x,{});
figure;
imagesc(y);
colorbar;


figure;
plot(L{2}.E,'r');

hold on;
plot(L{3}.E);


