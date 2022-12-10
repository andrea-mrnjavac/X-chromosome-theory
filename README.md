# SLiM simulation examples

## autosomal loci, beneficial mutations
```
initialize() {
        initializeMutationRate(1e-6);

        initializeMutationType("m1", h, "f", 0.03);

        initializeGenomicElementType("g1", m1, 1.0);

        initializeGenomicElement(g1, 0, 99);
        initializeRecombinationRate(0.5);

        initializeSex("A");
}

// create a population of 100 individuals
1 {
        sim.addSubpop("p1", 100);
}

// output 

1005000 late() { sim.outputFixedMutations("~/slim_results/auto_pos_"+h+"_INDEX_0.03s_100l_100i.out.fixed"); }
```

## hemizygous X loci

```
initialize() {
        initializeMutationRate(1e-6);

        initializeMutationType("m1", h, "f", 0.03);

        initializeGenomicElementType("g1", m1, 1.0);

        initializeGenomicElement(g1, 0, 99);
        initializeRecombinationRate(0.5);

        initializeSex("X",1);
}

// create a population of 100 individuals
1 {
        sim.addSubpop("p1", 100);
}



// output 

1005000 late() { sim.outputFixedMutations("~/slim_results/hemi_X_pos_"+h+"_INDEX_0.03s_100l_100i.out.fixed"); }
```

## diploid X loci

```
initialize() {
        initializeMutationRate(1e-6);

        initializeMutationType("m1", h, "f", 0.03);

        initializeGenomicElementType("g1", m1, 1.0);

        initializeGenomicElement(g1, 0, 99);
        initializeRecombinationRate(0.5);

        initializeSex("X",h);
}

// create a population of 100 individuals
1 {
        sim.addSubpop("p1", 100);
}


// output 

1005000 late() { sim.outputFixedMutations("~/slim_results/diplo_X_pos_"+h+"_INDEX_0.03s_100l_100i.out.fixed"); }
```
