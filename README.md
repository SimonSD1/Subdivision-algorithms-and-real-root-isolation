# CCA-project : Subdivision algorithms and real root isolation
Simon Sepiol-Duchemin, Joshua Setia

## Compilation and Execution    
To compile the project (from src folder) :  
```bash
make
```
- The above command will compile all .c files of the project, and provide one executable per files containing a main() function.    
- The executables will be located in the src folder.  

## File/Folders organization
```
.  
├── DATA  
│   ├── FixedCoeffSize_ChangingDegree  
│   └── FixedDegree_ChangingCoeffSize  
├── README.md  
├── Reports  
└── src  
    ├── EfficiencyTests  
    │   └── Results  
    ├── HeaderFiles  
    ├── Implementations  
    ├── ValidityTests  
    └── makefile
```

- main() functions are only in files located in folders EfficiencyTests and ValidityTests.  
- Results folder contains the result datas and graphs (in png).   
- Graphs are created with a python script located in folder EfficiencyTests, and should be generated directly in .c main() functions (not from terminal).
- DATA contains randomely generated polynomials that can be reused for future tests
