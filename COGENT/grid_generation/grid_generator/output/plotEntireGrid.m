files = dir('coords*');

for i=1:length(files)
    eval(['load ' files(i).name ' -ascii']);
    eval(['plot_coords(' files(i).name ')']);
end

print -dpng picture.png

