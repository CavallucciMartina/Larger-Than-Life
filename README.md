# Larger Than Life

## Istruzioni per compilazione

Con il comando `make` viene eseguito il file Make che automaticamente contiene
le regole di compilazione sia per OpenMp che CUDA.

## Istruzioni per esecuzione 

###Versione OpenMp
Lanciare, ad esempio, il seguente comando: `./omp-ltl 1 3 3 3 4 10 input.pbm output.pbm`
che corrisponde a ./omp-ltl (R =) 1 (B1 =) 3 (B2 =) 3 (D1 =) 3 (D2 = )4 (nsteps) 10 input_file output_file

###Versione CUD
Lanciare, ad esempio, il seguente comando: `./cuda-ltl 1 3 3 3 4 10 input.pbm output.pbm`
che corrisponde a ./cuda-ltl (R =) 1 (B1 =) 3 (B2 =) 3 (D1 =) 3 (D2 = )4 (nsteps) 10 input_file output_file
