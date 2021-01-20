# Bullseye
Bullseye plotter for converting interaction maps to Manhattan space

<img src="https://github.com/5centmike/Bullseye/blob/main/vc5C7.png" width="400" height="400">

This tool converts Euclidean distance plots to Manhattan or L1 Norm circular plots from the center point. The original use for this operation was visualizing point-to-point interaction among and between chromosomes in the nucleus. Conformation information in the genome provides a dense 2D matrix of interaction frequencies between all loci in the genome. In mammalian genomes strong points of interactions represent loops, generally formed between covergently oriented binding sites of the protein CTCF. These loops strongly influence the interactions among the adjacent loci along the looping polymers such that the effective distance between two loci is better described as a function of the sum of their distances to the loop.

![loops](https://github.com/5centmike/Bullseye/blob/main/loop.png)

This tool reshapes a 2D heatmap to produce a circular plot where the Euclidean distance within the circular plot is equal to the Manhatten distance from the center point in the original. In the figure below the data points in the blue diamond have been remapped into the blue circle.

![plots](https://github.com/5centmike/Bullseye/blob/main/plots.png)

This tool allows for a better visualization and understanding of the dynamics of polymer looping interactions, and I assume is likely applicable to other polymer interaction analyses beyond genomics.
