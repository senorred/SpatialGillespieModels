# SpatialGillespieModels
2D and 3D models of transcription/translation in confined/crowded areas

##### Crowded and Confined Spatial Gillespie Models 
These simulations take a standard gene expression model and add 2D or 3D migration of mobile components (here termed "ribos", but not neccessarily representative of translation). The intent here was to simulate localization of gene expression resources around a gene expression center, while modulating spatial access through confinement (reflective boundaries in different environment sizes) and crowding (obstacles in the for of occupied pixels or voxels).

By watching the distribution of gene expression resources over time, you can see distinct regions of influence form around stationary gene expression centers as resources are acquired from local areas. I called these regions of influence around gene expression centers "watersheds", labelled by a distinct color in the experiment space as the simulation progresses. One interesting finding was seeing how highly crowded regions were strongly limited in their resource gathering by the spatial limitation of obstacles or boundaries in the local space.

![50x50 pixel watersheds](https://github.com/senorred/SpatialGillespieModels/blob/master/trial12_50x50.png)

I also called these "hippo plots" or a "hippoplotamus", because I thought the idea of gene resources being captured by adversarial fixed expression centers in a shared space was similar to the game "hungry hungry hippos". Turns out that particular game is not as universal as I thought and was not a very good didactic metaphor, despite my best efforts to make it work.
