# FLAMES Library for Vitis HLS
Flexible Linear Algebra with Matrix-Empowered Synthesis (for Vitis HLS)

Developed by [Wuqiong Zhao](https://wqzhao.org) and other contributors,
from LEADS, Southeast University.

## Supported Versions
Since FLAMES is a modern library written using C++14 (with some C++17 features),
it only supports Vitis HLS 2020 or later where the CCC compiler version is 6.2.0.
Notably, Vivado HLS is not supported.

## Usage
FLAMES is a **header-only** library, so you can easily use it by first cloning it
under you project root with
```
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
