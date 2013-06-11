% get the histogram 
[count edges mid loc] = histcn([(simK'-kss)/kss*100 (simZ'-1)*100 ]);

% make a grid for plotting 
[X Y]=ndgrid(edges{1}, edges{2}); 
X=X(:); Y=Y(:);

scatter3(X,Y,count(:))

% calculate sizes so the most dense cell gets a value of 100 
% also convert from volume to "area" (as if drawing a sphere with 
% the right volume and cross-sectional area s) 
s_scale = 100/(max(count(:))^(2/3)); 
s = count(:).^(2/3) * s_scale; 
% convert any zeros to small numbers for scatter3 
s(s==0)=realmin;

% plot the densities 
fh=figure(); 
set(fh, 'Renderer', 'OpenGL'); % faster drawing 
scatter(X, Y, s, 'filled');