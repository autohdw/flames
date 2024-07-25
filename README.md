# FLAMES Library for Vitis HLS
Flexible Linear Algebra with Matrix-Empowered Synthesis (for Vitis HLS)

Developed by [Wuqiong Zhao](https://wqzhao.org) and other contributors,
from LEADS, Southeast University.

## Citation

If you find FLAMES useful, please cite our paper: [**Flexible High-Level Synthesis Library for Linear Transformations**](https://wqzhao.org/assets/zhao2023flexible.pdf), accepted by IEEE TCAS-II.
```bibtex
@article{zhao2024flexible,
  title     = {Flexible High-Level Synthesis Library for Linear Transformations},
  author    = {Zhao, Wuqiong and Li, Changhan and Ji, Zhenhao and Guo, Zhichen and Chen, Xuanbo and You, You and Huang, Yongming and You, Xiaohu and Zhang, Chuan},
  journal   = {{IEEE} Transactions on Circuits and Systems {II}: Express Briefs},
  volume    = {71},
  number    = {7},
  pages     = {3348--3352},
  year      = {2024},
  month     = jul,
  publisher = {IEEE}
}
```

## Supported Versions
Since FLAMES is a modern library written using C++14 (with some C++17 features),
it only supports Vitis HLS 2020 or later where the GCC compiler version is 6.2.0.
Notably, Vivado HLS is not supported.

## Usage
FLAMES is a **header-only** library, so you can easily use it by first cloning it
under you project root with
```shell
git clone https://github.com/autohdw/flames.git
```

### Core Modules
You can include the required header
```cpp
#include "flames/core.hpp"
```
All classes and functions are under the `flames` namespace, so you can use
```cpp
using namespace flames;
```
to access classes `Mat`, etc. directly instead of `flames::Mat`.

### Additional Insights
You can find more information about FLAMES in [`FLAMES_Insight.pdf`](https://flames.autohdw.com/FLAMES_Insight.pdf).

## License
The FLAMES is open source and distributed by an Apache license (v2.0).
Please [cite our paper](#citation) if you use FLAMES in your research.
