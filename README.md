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

### Artifact

HE_Mesh is hosted as an Maven artifact via Jitpack:

* Step 1. Add the JitPack repository to your build file

```
<repositories>
 <repository>
     <id>jitpack.io</id>
     <url>https://jitpack.io</url>
 </repository>
</repositories>
 ```
 
* Step 2. Add the dependency
 
 ```
<dependency>
    <groupId>com.github.micycle1</groupId>
    <artifactId>HE_Mesh</artifactId>
    <version>1.0.0</version>
</dependency>
 ```
