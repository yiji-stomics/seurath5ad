# seurath5ad

This package extends Seurat to enable the loading of raw H5AD files generated by other platforms or tools. Only the `obs` and spatial information from the H5AD files will be imported.

## Usage
```r
library(seurath5ad)

#Example usage  
#Load a raw H5AD file  
seurat_object <- LoadRawH5ad("path/to/your/file.h5ad")
```
## Installation

You can install the development version of seurath5ad from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("yiji-stomics/seurath5ad")
```

## Noted
When loading certain H5AD files, error messages may indicate that Seurat requires matrices in “double” format. It is essential to first convert the matrices to double format using Python.
eg.

```python
import numpy as np
data = data.X.astype(np.double)
```
