function setup ()

compinf = mex.getCompilerConfigurations ('C++') ;	% checking wether a C++ compiler configuration is available

fprintf ('Changing the current directory to ''../src'' ... ') ;

cd '../src' ;						% go to source directory

fprintf ('ok\nCompiling the pcaPP package ... ') ;
							% compiling the C++ sources
mex -DMATLAB_MEX_FILE -llibmwblas -llibmwlapack            ... 
   pcaPP.cpp L1Median_HoCr.cpp L1Median_VardiZhang.cpp     ...
   ML_meal.cpp ML_package.cpp ML_passrng.cpp outSDo.cpp    ...
   PCAgrid.cpp PCAproj.cpp qnn.cpp smat.cpp

fprintf ('ok\nCopying the ''pcaPP.mex*'' file(s) to the ''../matlab'' directory ... ') ;

copyfile ('pcaPP.mex*', '../matlab')			% copying the mex file(s) to the matlab directoy

fprintf ('ok\nChanging the current directory back to ''../matlab'' ... ')

cd '../matlab' ;					% changing back to the matlab directory

fprintf ('ok\n\n  Successfully compiled the pcaPP package for Matlab!\n') ;
