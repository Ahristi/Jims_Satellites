import csv
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon
import matplotlib.image as mpimg
import numpy as np
from Satellite import Satellite



def mappingHeatmap():
    error = []
    patches = []

    # CSV FORMAT : BoundingBox_corner_Long1, BoundingBox_corner_Lat1, BoundingBox_corner_Long2, BoundingBox_corner_Lat2, mapping_error
    with open('bounding/ALL_RESULTS.csv', mode='r') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            error.append(float(row[4]))
            # Create polygon based on long, lats
            x1 = float(row[0])
            y1 = float(row[1])
            x3 = float(row[2])
            y3 = float(row[3])

            x2 = x1
            y2 = y3
            x4 = x3
            y4 = y1

            polygon = Polygon([(x1,y1), (x2,y2), (x3,y3), (x4,y4)], closed = True)
            patches.append(polygon)

    error = np.array(error)

    fig, ax = plt.subplots(figsize=(10, 8))
    img = mpimg.imread('BlueMarble.png')
    ax.imshow(img, extent=[-180, 180, -90, 90])
    ax.grid(True)
    ax.set_xlabel('Longitude (deg)')
    ax.set_ylabel('Latitude (deg)')
    ax.set_xlim([135, 155])  # Longitude bounds for NSW
    ax.set_ylim([-38, -28])  # Latitude bounds for NSW
    ax.set_title('Mapping Error of Satellites')
    norm = plt.Normalize(0, 50)
    collection = PatchCollection(patches, cmap = 'plasma', norm = norm)
    collection.set_array(error)
    boxes = ax.add_collection(collection)
    fig.colorbar(boxes, label = 'Mapping Error [m]')
    ax.autoscale_view()
    plt.show()


if __name__ == "__main__":
    mappingHeatmap()
    