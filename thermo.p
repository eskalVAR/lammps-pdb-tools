set terminal dumb 90,30
datafile = "thermo_data"
firstrow = system('head -1 '.datafile)
set xlabel word(firstrow,1)
plot datafile using 1:5 title word(firstrow,5)

