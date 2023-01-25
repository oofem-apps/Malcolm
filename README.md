# MALCOLM
[![License: LGPL v2.1](https://img.shields.io/badge/License-LGPL%20v2.1-blue.svg)](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)


The aim of program MALCOLM is to assist structural engineers with the design of columns with multi-spiral reinforcment (MSR). This is done by means of an interaction diagram which presents an envelope of all admissible combinations of internal forces the cross section can sustain. Unfortunately, in currently valid design codes the MSR columns are not recognized which prohibits making use of the superb structural performance originating from passively confined concrete; in consequence the strength of MSR columns exceeds the resistance of standard columns by approx. 25%. Therefore, the approach implemented in MALCOLM should result into immediate material savings as the potential of the structural materials is utilized more efficiently.

In MALCOLM the algorithm is wrapped up in an accessible and friendly graphical user interface. Using this approach, key input parameters--the topology of the reinforcement and the material properties of concrete and steel--can be tuned to deliver the interaction diagram that satisfies all design load combinations. 
The resulting interaction diagram can be subsequently compared with the diagram constructed from nonlinear simulations, which is done here in the open source OOFEM finite element package. For more information, please refer to the [manual](malcolm_manual.pdf).

Program MALCOLM can redistributed it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.

For educational and research purposes, license for using MALCOLM is free of charge.
For commercial purposes, please contact program developer.

Copyright (C) 2022 Petr Havlasek


![malcolm_gif](malcolm_gif.gif)


### Pre-requisites

* In addition to common Python libraries it is essential to install PySide2 (GUI)
```
$ pip3 install pyside2
```
* OOFEM configuration with following flags turned ON:
```
USE_DSS (direct sparse solver)
USE_OPENMP (parallel computation using multiple threads)
USE_PYBIND_BINDINGS (python bindings for external calling OOFEM from MALCOLM)
```
## Running MALCOLM

```
$ python3 malcolm.py
```

## Documentation
MALCOLM documentation/user manual is available [here](malcolm_manual.pdf).


## Authors
[Petr Havlásek](mailto:petr.havlasek@cvut.cz), [ORCID: 0000-0002-7128-3664](https://orcid.org/0000-0002-7128-3664)<br/>
Alena Plačková (MALCOLM logo)

## Acknowledgments
Financial support for this work was provided by the Technology Agency of the Czech Republic (TAČR), project number TM01000059 (Reducing material demands and enhancing structural capacity of multi-spiral reinforced concrete columns - advanced simulation and experimental validation).