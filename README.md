# haplot
An easy-to-use python package can be used for omics visual analysis

## Installation
```bash
pip install haplot
```

## Usage
```python
import haplot as hp
```

## Functions
### 1. plot_heatmap
```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from haplot.chart import LDHeatmapPlot

np.random.seed(123)
# generate random R^2 value for ld heatmap plot
df = pd.DataFrame(
    data=np.random.uniform(0, 1, (10, 10)),
    index=['rs{}'.format(i + 1) for i in range(10)],
    columns=['rs{}'.format(i + 1) for i in range(10)]
)

fig, ax = plt.subplots(figsize=(5, 5))
LDHeatmapPlot(df)
plt.show()
```

![test5](tests/test5.png)
