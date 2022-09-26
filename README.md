# Rigid body physics
Code by Charles Gannon cgannon@ucmerced.edu

The author (s) make no guarantees to accuracy or precision of this code.

The author (s) do not claim ownership of nor claim to have contributed to any dependencies.

## Description
A rigid body physics engine designed for educational purposes.

## Dependencies

### Eigen

This code uses the 3rd party Eigen c++ library.

https://eigen.tuxfamily.org/index.php?title=Main_Page

This Eigen library is used under the MPL 2.0 license.

https://www.mozilla.org/en-US/MPL/2.0/

### Autodiff
This code uses the 3rd party Autodiff c++ library.

https://autodiff.github.io/

This Autodiff library is used under the MIT liscense.

https://www.mit.edu/~amini/LICENSE.md

## Status

Currently this project is a work in progress and is currently (very) incomplete.

The project is likely riddled with bugs and glitches and inaccuracies.

Therefore, using this project for anything but example code is not recommended (at least for now).

## Philosophy

### General Philosophy
Physical simulations were historically one of the drivers of the improvement in software and computers in general.
As such physics simulation software has been done before and has been done much better than one person could hope to accomplish.
Therefore, this software will not aim to be more numerically precise, faster or feature rich than other simulation software. 
Instead, the aim of this project is to create physics software for education in two senses. Firstly, the software
will aim to provide intuitive physics simulation for education, and will target it's feature set towards education. Secondly the source 
code and documentation of the software will aim to be educational itself.

### Design Philosophy
- Beauty in simplicity, complexity only when necessary
- Code should imply its own functionality, while the documentation makes functionality explicit
- Readability and accuracy first, performance second
- State should be separate from logic
- Reinventing the wheel is not necessary
- Write once use many times

## Goals

### Source Code
- The source code is readable and reusable by a 3rd party
- The source code is sufficiently documented such that a 3rd party would be able to understand the function of code
- Any special cases are clearly noted in comments and documentation
- The source code maintains a uniform format (IE bracket placement, tab spacing, ...)
- The source code adheres to the philosophy outlined in the philosophy section of this document
- The source code maintains a high degree of reusability and extensibility

### Physics
- All critical physics / mathematical equations are documented with references to external sources.
- All physics derivations are accurate.
- Various test are set up such that simulations are compared against analytical and numerical results.

### Miscellaneous
- Simulations are reasonably performant even for the most demanding simulations a user should reasonably expect
- Code should be committed and pushed to git at the end of every workday
- All documentation should be provided in nonproprietary formats

## Outline

The outline for this project. Phase 1 has detailed tasks. Later phases only have vague goals.

### Phase 1 = Physics core: part 1 (Point particles)
Physics of point particles
- [x] Write core class data structures
- [x] Implement Euler method for 2nd order ODEs
- [x] Write core evolution code
- [ ] Document everything up to this point
- [x] Implement logger
- [x] Write state recorder / basic I/O code
- [x] Implement constraint forces
- [z] Implement Runge-Kutta method for 2nd order ODEs
- [ ] Document everything up to this point
- [ ] Write tests

### Phase 2 - Physics core: part 2(Bodies)
Physics of non-point bodies
- Research implementing non-point bodies 
- Implement rotation
- Implement collisions
- Implement tests

## Phase 3 - Refactoring
- Compile physics core to library, separate executable
- Implement parallelization
- Optimize previously written code

## Phase 4 - Metrics
Implement various metrics intended for education like tracking kinetic / potential energy.

## Phase 5 - Command line interface
Design and implement command line assembly for interfacing with physics engine as well as 
implement file I/O system.

## Phase 6 - Web / UI
Compile to physics core to web assembly. Implement graphical user interface through browser.

## LISENCE

This software is provided under the MIT liscense.
https://www.mit.edu/~amini/LICENSE.md

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.