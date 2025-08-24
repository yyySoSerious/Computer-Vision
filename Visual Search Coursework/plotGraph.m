figure;

%global color histogram
rgb_qlevel = [2, 3, 4, 5, 6, 7, 8, 10, 12];
rgb_map = [0.173, 0.193, 0.210, 0.214, 0.210, 0.212, 0.208, 0.202, 0.197];

%spatial grid
sgrid_qdivs = [2, 3, 4, 5, 6, 7, 8, 10, 12, 14];
sgrid_map = [0.186, 0.205, 0.210, 0.211, 0.209, 0.209, 0.208, 0.208, 0.205, 0.204];

%eoh
eoh_qintervals = [2, 4, 8, 16, 17, 18, 19, 20, 32, 64];
eoh_map = [0.143, 0.154, 0.213, 0.226, 0.226, 0.227, 0.225, 0.225, 0.219, 0.213];

%eoh with colour
eohWithColor_qintervals = [2, 4, 8, 16, 17, 18, 19, 20, 32, 64];
eohWithColor_map = [0.210, 0.213, 0.233, 0.236, 0.242, 0.240, 0.239, 0.235, 0.232, 0.228];

x = eohWithColor_qintervals;
y = eohWithColor_map; 
p = plot(x, y);
title('MAP values for EOH with Color ');
xlabel('{$q$} angular intervals', 'interpreter','latex', ...
    'fontsize', 12);
ylabel('MAP', 'interpreter','latex', 'fontsize', 12);