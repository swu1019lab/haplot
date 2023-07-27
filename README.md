# haplot
![Static Badge](https://img.shields.io/badge/github-haplot-blue)

An easy-to-use python package can be used for omics visual analysis

## Dependencies
- python >= 3.10
- numpy >= 1.22.3
- pandas >= 1.4.2
- matplotlib >= 3.6.2
- scipy >= 1.9.3
- geopandas >= 0.13.2


## Usage
```python
import haplot as hp
```

## Functions
### 1. LDHeatmapPlot
```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from haplot.chart import LDHeatmapPlot

np.random.seed(123)
# generate random R^2 value for ld heatmap plot
# snp id and snp position
arrays = [['rs{}'.format(i + 1) for i in range(10)],
          np.sort(np.random.randint(1, 1000000, 10))]
index = pd.MultiIndex.from_arrays(arrays, names=('rs_id', 'pos'))

df = pd.DataFrame(
    data=np.random.uniform(0, 1, (10, 10)),
    index=index,
    columns=index
)

fig, ax = plt.subplots(figsize=(5, 5))
LDHeatmapPlot(df, plot_value=True, plot_snp=True, plot_diag=False, cmap='Reds')
plt.show()
```

![test5](tests/test5.png)
