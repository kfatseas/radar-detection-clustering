# radar-detection-clustering

Provides a minimal class to load data from the recordings that have been created with our radar-camera setup. The Recording class reads the radar settings and handles all the digital signal processing required to extract detections from each radar frame.

To demonstrate the clustering of radar detections, we utilize the DBSCAN algorithm. The results are visualized in both the range-doppler map and the cartesian plane after transforming the polar coordinates to cartesian. Keep in mind that the detections which do not belong to any cluster (noise) do not appear on the plots. 

### Todo

- [ ] Try transforming the radar detection coordinates to cartesian before clustering
- [ ] Improve the range-doppler map filtering, especially for the static objects
