

This repository tracks work-in-progress research into a model of antigenic drift in Sars-CoV-2 by Peter Jentsch, Finlay Maguire, and Samira Mubareka.

Please email `peter_jentsch@sfu.ca` if you have any questions.

## To run
Ensure you have the latest stable version of Julia installed and in your `PATH`. This code tested on Linux only.

```
git clone https://github.com/pcjentsch/genomic_landscape_evolution.git
cd genomic_landscape_evolution
julia --project=SarsEvoModel
```

then

```
]instantiate
using SarsEvoModel
main() 
```


