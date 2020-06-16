%Loading file 'InputData.csv" to Workspace

dane_y = dlmread('InputData.csv', ',', [2,1,8,30]);
dane_x1 = dlmread('InputData.csv', ',', [12,1,18,30]);
dane_x2 = dlmread('InputData.csv', ',', [22,1,28,30]);

y = dane_y(:);
[t, n] = size(dane_y);
x = [dane_x1(:) dane_x2(:)]; 