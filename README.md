# FLAMES Library for Vitis HLS
Flexible Linear Algebra with Matrix-Empowered Synthesis (for Vitis HLS)

Developed by [Wuqiong Zhao](https://wqzhao.org) and other contributors,
from LEADS, Southeast University.
This is part of the [AutoHDW](https://autohdw.com) project.

## Key Features
- Header-only Library *(Super Easy to Install and Use)*
- Templated-based C++ Design:
  - `Mat` and `Tensor` Classes
  - `MatRef` and `MatView` for Reference and View
  - (Relatively) Modern C++ Implementation (C++14/17)
- Hardware-friendly Design:
  - Optimized RAM Usage
  - Configurable Parallelism
  - Optimized Matrix Operations
- Fully Open Source

Read [our paper](#citation) for details ([open access PDF](https://wqzhao.org/assets/zhao2024flexible.pdf)).

## Example
Linear algebra algorithm implementation made easy by FLAMES!
Here is the NSA algorithm and the corresponding implementation using FLAMES:

| Step | Algorithm                                    | FLAMES C++ Implementation          |
| :--: | -------------------------------------------- | ---------------------------------- |
| 1    | $\mathbf{D} = \mathbf{A} \circ \mathbf{I}$   | `auto D = A.diagMat_();`         |
| 2    | $\mathbf{E} = \mathbf{A} - \mathbf{D}$       | `auto E = A.offDiag_();`         |
| 3    | $\mathbf{D}_I = \mathbf{D}^{-1}$             | `auto D_I = D.inv();`              |
| 4    | $\mathbf{P} = -\mathbf{D}_I \mathbf{E}$      | `auto P = -D_I * E;`               |
| 5    | $\mathbf{X} = \mathbf{P}$ (Iter. 1)          | `auto X = P_ = P;`                 |
| 6    | $\text{for } i = 2, \dots, n$                | `for (int i = 2; i <= n; ++i) {`   |
| 7    | $\quad \mathbf{P}^i = \mathbf{P}^{i-1} \mathbf{P}$ | `    P_ *= P;`               |
| 8    | $\quad \mathbf{X} = \mathbf{X} + \mathbf{P}^i$     | `    X += P_;`               |
| 9    | $\text{end}$                                 | `}`                                |
| 10   | $\mathbf{A}^{-1} = \mathbf{X} \mathbf{D}_I + \mathbf{D}_I$ | `A_inv = X * D_I + D_I;` |

*From Table 2 in [our paper](#citation), code available at [`examples/mat-inv-nsa`](examples/mat-inv-nsa).*

## Citation

If you find FLAMES useful, please cite our paper: [**Flexible High-Level Synthesis Library for Linear Transformations**](https://ieeexplore.ieee.org/document/10437992), *IEEE TCAS-II*, vol. 71, no. 7, pp. 3348-3352, Jul. 2024. [[IEEE Xplore](https://ieeexplore.ieee.org/document/10437992)] [[PDF](https://wqzhao.org/assets/zhao2024flexible.pdf)] [[DOI](https://doi.org/10.1109/TCSII.2024.3366282)]
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
