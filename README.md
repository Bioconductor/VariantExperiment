# VariantExperiment

Herve: 
You need to stick the assay data stored in the GDS file into the SE. For this you need a GDS back-end for DelayedArray. I can provide more details if needed but you'll basically need to implement 2 classes for this e.g. GDSArraySeed and GDSArray. Note that this is a project on its own and should preferably go in a dedicated package (e.g. GDSArray).
