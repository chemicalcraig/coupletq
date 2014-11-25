#Exciton transfer calculation

calculation
  type pert
  ewindow 0.5
  molecules 3
end

molecules
  states 2
  charges 
    1 0 file.dat
    0 1 file.dat
    1 1 file.dat
  end
  move
    x -3. -3. 1
  end
  rotate
    x M_PI/2
  end
  tddft file.out
end

dynamics
  start a
  finish b
  steps c
  increment d
  initial
    0 1 1
  end
  output
    populations (0,0) (1,1)
    file file.dat
  end
end

