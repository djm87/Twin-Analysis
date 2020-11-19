i=2
j=2
baseDir=pwd;
temps={'625C','800C'}
direcitons={'ND','TD','RD'}
strain_lvls={'initial','5pct strain','20pct strain'}
pipeName='run_pipe'
pipePath=fullfile(baseDir,temps{i},direcitons{j},strain_lvls{1},pipeName)
run(pipePath);
