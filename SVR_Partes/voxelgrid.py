import numpy as np
import pandas as pd
from pyntcloud import PyntCloud

positions = np.random.rand(1000, 3) * 10
positions -= positions.mean(0)

points = pd.DataFrame(
    positions.astype(np.float32), 
    columns=['x', 'y', 'z'])

points["red"] = (np.random.rand(1000) * 255).astype(np.uint8)
points["green"] = (np.random.rand(1000) * 255).astype(np.uint8)
points["blue"] = (np.random.rand(1000) * 255).astype(np.uint8)

cloud = PyntCloud(points)
lines = [
    {
        "color": "red",
        "vertices": [[0, 0, 0], [10, 0, 0]]
    },
    {
        "color": "green",
        "vertices": [[0, 0, 0], [0, 10, 0]]
    },
    {
        "color": "blue",
        "vertices": [[0, 0, 0], [0, 0, 10]]
    }
]

cloud.plot(polylines=lines)