# Goal
* the main goal of this repo is to identify genes in the 3 dpf zebrafish opercle that are spatially informative, providing both potential biological meaning (genes which are restricted to one area might be involved in local functions) and spatail meaning that can inform better integration of single cell and spatial (atlas) data using Tangram or similar tools.

# Context
* using image registration, created a 3D gene expression atlas of the zebrafish opercle where about 30 genes were mapped into a common coordinate space.
* I used tangram to impute expression of unmapped genes into atlas space.
    * I used zscape single cell data to do this. Specifically, cells from the 3 dpf osteoblast cluster.
# Problem
* The main problem I have with Tangram predictions is that they are really most powerful for where we already have data.
* We can see this in the leave on out tests where it does a good job for genes that are within the core area or close to it, but a much poorer job for regions peripheral to it.
* This is problematic on both sides. (single cell and atlas data)
    * The single cell data might not contain as many of these cells.
	* Our atlas data also contains little information in these areas.
        * most genes map to a central osteoblastic region where committed osteoblastic marker genes are expressed.
        * This distribution seems logical since atlas genes were initially selected based on OB marker genes.
* There's not much for getting around the fact that we don't have much spatial information at this point for the surrounding area. 

# Plan (conceptual)
* We have an expectation that peripheral cells are part of a neural crest derrived mesenchymal condensation that eventually gives rise to the opercle. 
* Can we use this information to recluster the singel cell data to capture the whole condensation at 3 dpf.
* challenges I see are: 
	how are we going to separate out the cells in this condensation from any other condensataion?
	What are the markers of this condensation?

# Plan (potential implementation)
Maybe one way around this is to use trajectory analysis to figure out the progenitor populations of the annotated osteoblasts and use that to expand the population of cells we try to find markers from.
From this, we can use tools like spapros to find coarse markers (cell type and variation markers) in the broader OB+progenitor population.
We can then map these genes onto our atlas to inform better tangram predictions of spatially restricted genes. 

