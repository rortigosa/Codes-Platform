function out = preplot_path()
 paths = {'c:\Program Files\Tecplot\Tec360 2011R2\bin\preplot.exe', ...
     'c:\Program Files\Tecplot\Tec360 2009\bin\preplot.exe', ...
     '/usr/bin/preplot.linux64.26',...
     '/usr/local/tecplot/bin/preplot',...
     };
 out = '';
 for s = paths
 if exist(s{1},'file')
     out = s{1};
 end
 end
 