<img src="https://user-images.githubusercontent.com/21171362/232594749-6a6beacc-6c59-4e40-ae62-e7026f355742.png" align="right" width="200"/>

# cellula

**`cellula`** is a simple R-based pipeline for single cell RNA-seq processing with a number of methods for integration and identity assignment.

You can read all about it on the companion [website](https://gdagstn.github.io/cellulaweb)!

## Changelog

-   13/3/2025:\
	Added heatmap plotting and proportion plotting functions

-   8/3/2025:\
	Added pseudobulk DE functions for analysis and plotting 
	
-   4/3/2025:\
    Added integration diagnostics

-   29/10/2024:\
    Added stacked violin plots and multipanel dimensionality reduction plots

-   29/1/2024:\
    Added additional methods for integration: scMerge2 from `scMerge` and STACAS from `STACAS`
    
-   30/7/2023:\
    Added a function to build references for the Jaitin method from a SCE. Fixed minor bugs in identity assignment

-   23/7/2023:\
    Added several options for color palettes by overhauling color management.

-   21/7/2023:\
    Moved many optional dependencies in Suggests in DESCRIPTION, added cell cycle estimation and visualization

-   15/6/2023:\
    Added a new fast method for identity assignment based on a reference (Jaitin method)

-   3/6/2023:\
    Added interface to monocle3 trajectory inference and a new UMAP embedding
