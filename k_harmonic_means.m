clc;
clear;
close all;

% K-Harmonic Means
% Author: Tanadon Ratsameephen

load fcmdata.dat

k = 2;
c = 0.01;
dataSize = size(fcmdata, 1);
randomcentersIdx = randperm(dataSize, k);
centers = fcmdata(randomcentersIdx, :);

p = 2;
maxIter = 5;

khm = 0; % khm
for i = 1:dataSize
    denominator = 0;
    for j = 1:k
        denominator = denominator + (1 ./ sqrt(sum((fcmdata(i,:)-centers(j,:) + c).^2, 2)).^p);
    end
    khm = khm + (k./denominator);
end

khmNew = khm;

while true
    m = zeros(k, dataSize);
    for i = 1:dataSize
        denominator = 0;
        for j = 1:k
            denominator = denominator + sqrt(sum((fcmdata(i,:)-centers(j,:) + c).^2, 2).^(-p-2));
        end

        for j = 1:k
            % Numerator fixed.
            numerator = sqrt(sum((fcmdata(i,:)-centers(j,:) + c).^2, 2).^(-p-2));
            m(j,i) = numerator / denominator;
        end
    end

    w = zeros(1, dataSize);
    for i = 1:dataSize
        denominator = 0;
        for j = 1:k
            denominator = denominator + sqrt(sum((fcmdata(i,:)-centers(j,:) + c).^2, 2).^(-p-2));
        end

        w(i) = denominator / denominator.^2;
    end
    
    %%
    khmOld = khmNew;
    
    khm = 0; % khm
    for i = 1:dataSize
        denominator = 0;
        for j = 1:k
            denominator = denominator + (1 ./ sqrt(sum((fcmdata(i,:)-centers(j,:) + c).^2, 2)).^p);
        end
        khm = khm + (k./denominator);
    end
    
    khmNew = khm;
    
    if(khmNew/khmOld < 0.6)
        break;
    end
    %%

    for j = 1:k
        numerator = zeros(size(fcmdata(i,:)));
        denominator = 0;
        for i = 1:dataSize
            numerator = numerator + m(j, i) * w(i) .* fcmdata(i,:);
            denominator = denominator +  m(j, i) * w(i);
        end

        centers(j,:) = numerator ./ denominator;
    end
    
end

maxU = max(m);
index1 = find(m(1,:) == maxU);
index2 = find(m(2,:) == maxU);

plot(fcmdata(index1,1),fcmdata(index1,2),'ob')
hold on
plot(fcmdata(index2,1),fcmdata(index2,2),'or')
plot(centers(1,1),centers(1,2),'xb','MarkerSize',15,'LineWidth',3)
plot(centers(2,1),centers(2,2),'xr','MarkerSize',15,'LineWidth',3)
hold off


