# radar-detection-clustering

Provides a minimal class to load data from the recordings that have been created with our radar-camera setup. The Recording class reads the radar settings and handles all the DSP processing that is required to extract detection from each radar frame.

To demostrate the clustering of radar detections, the DBSCAN algorithm is utilized. The results are visualized in both the range-doppler map but also the cartesian plane after tranforming the polar coordinates to cartesian. Keep in mind that the detections do not belong to any cluster (noise) are not shown. 


### Todo

- [ ] Try tranforming the radar detection coordinates to cartesian before clustering
