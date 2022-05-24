# -*- coding: utf-8 -*-
"""
Created on Thu May 12 23:00:33 2022

@author: Konstantinos
"""

import matplotlib.pyplot as plt
from PIL import Image
from scipy import signal
import numpy as np
import pandas as pd
import os


def mag(signal):
    return 20 * np.log10(abs(signal))

# Minimal recording class to load data from recording folder
class Recording():
    def __init__(self, directory):
        self.directory = directory
        self.settings, self.profiles, self.meta = self.parse_info(directory)
        self.noise_floor = self._noise_floor()

    def __len__(self):
        return self.meta['frames saved'].item()

    # Read radar settings and meta information
    def parse_info(self, directory):
        filename = os.path.join(self.directory, "Settings.csv")
        settings = pd.read_csv(filename, nrows=1)
        profiles = pd.read_csv(filename, skiprows=2)
        filename = os.path.join(self.directory, "meta.csv")
        meta = pd.read_csv(filename)
        return settings, profiles, meta

    # Returns the image taken from the camera of the setup
    def image(self, idx):
        camera_filename = os.path.join(self.directory, f'camera0 {idx}.jpg')
        try:
            image = Image.open(camera_filename)
        except FileNotFoundError:
            tmp = np.zeros((720, 1280, 3), dtype=np.uint8)
            image = Image.fromarray(tmp)
        return image

    # Read the file with the raw ADC values
    # Reashape the array so that it has the following shape
    # (virual_antennas, chirps, samples per chirp)
    def raw(self, idx):
        chirps = self.settings['chirps/tx'].item()
        samples = int(self.settings['all samples'].item() / self.profiles.iloc[0]['deci'])
        ant = self.settings['virtual antennas'].item()
        tx = self.settings['tx antennas'].item()
        rx = self.settings['rx antennas'].item()
        filename = os.path.join(self.directory, "raw " + str(idx) + ".dat")
        try:
            frame = np.fromfile(filename, dtype=np.uint16)
            frame = (frame / 2048.0) - 1.0
            frame = frame.reshape(rx, chirps, tx, samples)
            cube = np.empty((ant, chirps, samples), dtype=np.float64)
            for t in range(tx):
                for r in range(rx):
                    cube[t*rx+r] = frame[r, :, t]
        # If file is not there, return an array of ones with the proper shape
        except FileNotFoundError:
            cube = np.ones((ant, chirps, samples), dtype=np.float64)
        return cube

    # Returns the radar cube for the specified frame index after performing a 2D FFT.
    # The 2 FFTs are computing the range and velocity
    # Window is applied before each FFT to reduce spectrum leakage
    def rd_cube(self, idx):
        # Get Radar Cube from Raw File
        raw_cube = self.raw(idx)
        antennas, chirps, samples = raw_cube.shape
        # Compute Range FFT
        window = signal.windows.hann(samples)[np.newaxis, np.newaxis, :]
        raw_cube = raw_cube * window
        range_cube = np.fft.rfft(raw_cube, axis=2)
        range_cube /= samples * (1 / np.mean(window))
        range_cube = range_cube[:, :, 0:samples//2]
        # Compute Velocity FFT
        window = signal.windows.hann(chirps)[np.newaxis, :, np.newaxis]
        range_cube = range_cube * window
        range_doppler_cube = np.fft.fft(range_cube, axis=1)
        range_doppler_cube = np.fft.fftshift(range_doppler_cube, axes=1)
        range_doppler_cube /= chirps * (1 / np.mean(window))
        return range_doppler_cube

    # To find the noise floor and use it in order to threshold the range-doppler map
    # we compute an average of some (50 in this case) range-doppler maps.
    # Then we remove the static objects and compute the threshold value per column
    def _noise_floor(self):
        chirps = self.settings['chirps/tx'].item()
        rds = [mag(self.rd_cube(i)).max(axis=0) for i in range(50)]
        rds = np.array(rds)
        threshold = rds.mean(axis=0)
        i = int(chirps/2)
        threshold = np.delete(threshold, (i-2, i-1, i, i+1, i+2), axis=0)
        threshold = np.mean(threshold, axis=0)
        noise_floor = np.repeat(threshold[np.newaxis, :], chirps, axis=0)
        return noise_floor

    # Returns the range-doppler map which is computed by keeping the max value
    # of the 2DFFT processed cube along the antenna axis.
    # The shape of the output is (chirps, samples)
    def range_doppler(self, idx):
        rd_cube = self.rd_cube(idx)
        range_doppler = rd_cube.max(axis=0)
        return mag(range_doppler)

    # Returns a tuple with 
    # To filter the range doppler maps we set the threshold
    # at 10dB above the noise floor.
    # All cells of the range-doppler map with a higher value
    # are the detections that we gonna compute the DoA for
    def detections(self, idx):
        rd = self.range_doppler(idx)
        ys, xs = np.where(rd > rec.noise_floor + 10)  # The threshold is 10dB above the noise floor
        # For each detection perform the 3rd FFT on the antenna axis to find the direction of arrival(DoA)
        rd_cube = self.rd_cube(idx)
        antennas, chirps, samples = rd_cube.shape
        snapshots = rd_cube[:, ys, xs].T
        fft_size = 181
        window = signal.windows.chebwin(antennas, at=80)
        windowed_snapshots = snapshots * window[np.newaxis, :]
        padding = int((fft_size-antennas) / 2)
        padded_snapshots = np.pad(windowed_snapshots, [(0, 0), (padding, padding)], mode='constant')
        doa_spectrum = np.fft.fft(padded_snapshots, axis=1) / (antennas * (1 / np.mean(window)))
        doa_spectrum = abs(doa_spectrum)
        degs = np.argmax(doa_spectrum, axis=1)
        # DoA correction
        degs = (degs - 90) / 90
        degs = np.degrees(np.arcsin(degs))
        degs += 90
        # degs = 181 - degs
        
        # X axis index to range
        max_range = rec.settings['max distance'].item()
        samples = rec.settings['all samples'].item() // 2
        ranges = xs * (max_range / samples)
        
        # Y axis index to radial velocity
        max_velocity = rec.settings['max velocity'].item()
        chirps = rec.settings['chirps/tx'].item()
        velocities = ys * (max_velocity * 2) / chirps - max_velocity
        velocities *= -1
        
        # Merge the data we computed so far for each detection
        # points = np.stack((velocities, ranges, rd[ys, xs], degs)).T
        detections = np.stack((velocities, ranges, degs)).T
        return detections




######## Main
######
####
###
##
#

# Load a recording
rec_folder = 'C:/DataSet/bastilleveld_south_east_1'
rec = Recording(rec_folder)

# Choose a random index
idx = np.random.randint(0, len(rec))
# Or select a specific frame
# idx = 1899

# Show the image from the camera
rec.image(idx).show()

# Show unfiltered range-doppler map
rd = rec.range_doppler(idx)
plt.figure()
plt.imshow(rd)

# Get the detections as tuples of (velocity, range, DoA)
points = rec.detections(idx)


#### Clustering
###
##
#

from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler

# If values are scaled, DBSCAN params need to be adjusted
X = StandardScaler().fit_transform(points)
# X = points

db = DBSCAN(eps=.1, min_samples=5).fit(X)

# Get core samples for each cluster
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_

n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_ = list(labels).count(-1)

# Black removed and is used for noise instead.
unique_labels = set(labels)
colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, n_clusters_)]

max_distance = rec.settings['max distance'].item()
max_velocity = rec.settings['max velocity'].item()

# Plot clustered detections in Range vs Doppler plot
plt.figure(figsize=(10,5))
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = [0, 0, 0, 1]

    class_member_mask = labels == k

    xy = points[class_member_mask & core_samples_mask]
    plt.plot(
        xy[:, 1],
        xy[:, 0],
        "o",
        markerfacecolor=tuple(col),
        markeredgecolor="k",
        markersize=14,
    )

    xy = points[class_member_mask & ~core_samples_mask]
    plt.plot(
        xy[:, 1],
        xy[:, 0],
        "o",
        markerfacecolor=tuple(col),
        markeredgecolor="k",
        markersize=6,
    )
plt.axis((0, max_distance, -max_velocity, max_velocity))
plt.title("Range-Doppler - Estimated number of clusters: %d" % n_clusters_)
# plt.savefig('rd_clusters.png', dpi=200)
plt.show()


# Plot clustered detections on the cartesian plane
plt.figure(figsize=(10,5))
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = [0, 0, 0, 1]

    class_member_mask = labels == k

    xy = points[class_member_mask & core_samples_mask]
    # Polar to cartesian
    x = xy[:, 1] * np.cos(np.radians(xy[:, 2]))
    y = xy[:, 1] * np.sin(np.radians(xy[:, 2]))
    plt.plot(
        x,
        y,
        "o",
        markerfacecolor=tuple(col),
        markeredgecolor="k",
        markersize=14,
    )
    # Polar to cartesian
    xy = points[class_member_mask & ~core_samples_mask]
    x = xy[:, 1] * np.cos(np.radians(xy[:, 2]))
    y = xy[:, 1] * np.sin(np.radians(xy[:, 2]))
    plt.plot(
        x,
        y,
        "o",
        markerfacecolor=tuple(col),
        markeredgecolor="k",
        markersize=6,
    )
plt.axis((-max_distance, max_distance, 0, max_distance))
plt.title("Cartesian coordinates - Estimated number of clusters: %d" % n_clusters_)
# plt.savefig('ra_clusters.png', dpi=200)
plt.show()