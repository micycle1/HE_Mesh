[![](https://jitpack.io/v/micycle1/HE_Mesh.svg)](https://jitpack.io/#micycle1/HE_Mesh)
 [![Lines of Code](https://sonarcloud.io/api/project_badges/measure?project=micycle1_HE_Mesh&metric=ncloc)](https://sonarcloud.io/summary/new_code?id=micycle1_HE_Mesh) [![Code Smells](https://sonarcloud.io/api/project_badges/measure?project=micycle1_HE_Mesh&metric=code_smells)](https://sonarcloud.io/summary/new_code?id=micycle1_HE_Mesh)

# HE_Mesh

HE_Mesh, a Java library for creating and manipulating polygonal meshes. Aimed primarily at [Processing](http://processing.org/).

## Maven Fork
This fork hosts _HE_Mesh_ as a fully mavenised dependency, making it much easier to use the library in Maven Java projects or make source code modifications.

**NOTE**  this fork has not been extensively tested.

### Fork Details

* Where applicable, code that once wrapped other libraries now references the original libraries directly.
(e.g. straight-skeleton functionality now references camp-skeleton (`org.twak.camp`) directly rather than the `wblut.external.straightskeleton` wrapper).
* HE_Mesh code from `hemesh-external.jar` is now included in the source files (classes such as `WB_JTS`, `WB_QuickHull3D`, etc.).
