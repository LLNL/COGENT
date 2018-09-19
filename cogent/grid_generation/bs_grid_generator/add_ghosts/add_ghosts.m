function add_ghosts(in_file, out_file)

  blocks = {'lcore' 'rcore' 'lcsol' 'rcsol' 'lsol' 'rsol' 'lpf' 'rpf' 'mcore' 'mcsol'};

  fid = fopen(out_file,'w');

  for n = 1:length(blocks)
    [R, Z, nr, nr_extend, np, np_extend] = add_block_ghosts(in_file, char(blocks(n)));

    fprintf(fid,'%s %d %d %d %d\n', char(blocks(n)), nr, nr_extend, np, np_extend);
    for j = 1:np + 2*np_extend
      for i = 1:nr + 2*nr_extend
        fprintf(fid,'%20.13e %20.13e %20.13e %20.13e\n', R(i,j), Z(i,j), 0., 0.);
      end
    end
  end

  fclose(fid);

  return;
