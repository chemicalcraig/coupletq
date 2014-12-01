#Exciton transfer calculation

calculation
  type pert
  ewindow 0.5
  molecules 3
end

molecules
  states 2
  charges 
    1 0 densoutput-s10.dat
    0 0 densoutput-s0.dat
    1 1 densoutput-s1.dat
  end
  move
    x -3. -3. 1
  end
  rotate
    x M_PI/2
  end
  tddft pmimeth.out
end

molecules
  states 2
  charges 
    1 0 densoutput-s10.dat
    0 0 densoutput-s0.dat
    1 1 densoutput-s1.dat
  end
  move
    x -3. -3. 1
  end
  rotate
    x M_PI/2
  end
  tddft pmimeth.out
end

dynamics
  start a
  finish b
  steps c
  increment d
  init
    0 1 1
  end
  output
    populations (0,1) (2,3)
    file file.dat
  end
end

