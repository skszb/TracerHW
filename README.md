# Ray Tracing Homework

***

## Overview:
This is the final homework of _**Computer Graphics and Animation**_.

In this homework, based on the ray tracer framework provided by Professor O'Brien, 
I implemented:

1. ray-surface intersection (triangle surfaces and spheres)

2. Phong shading

3. Shadow

4. Reflection

5. Acceleration data structure (bounding volume hierarchy, 
using axis-aligned bounding boxes)

---

### Usage
The ray tracer takes 3 arguments, 2 required and 1 optional.

Required: a file of Neutral File Format as an input, outputs a file of Portable Pixmap Format. 

(e.g. trace [opts] input.nff output.ppm)

More information about Neutral File Format can be found at http://www.realtimerendering.com/resources/SPD/NFF.TXT

***

## Outputs:
### Benchmarks
I used the **teapot** and the **balls** as my test cases of shading, shadow, and reflection.

These files can be found at http://www.realtimerendering.com/resources/SPD/

<img src="OutputFiles\_benchmark_teapot.png" width=45%></img>  
<img src="OutputFiles\_benchmark_balls.png" width=45%></img>

### My Outputs
#### Basic shading

<img src="OutputFiles\teapot_shading.png" width=45%></img>
<img src="OutputFiles\balls_shading.png" width=45%></img>

#### Shading + shadow

<img src="OutputFiles\teapot_shading_shadow.png" width="45%"/></img>
<img src="OutputFiles\balls_shading_shadow.png" width="45%"/></img>

#### Shading + shadow + reflection

<img src="OutputFiles\teapot_shading_shadow_reflection.png" width="45%"/></img>
<img src="OutputFiles\balls_shading_shadow_reflection.png" width="45%"/></img>

#### Shading + shadow + reflection (more rays, num=5)

<img src="OutputFiles\teapot_shading_shadow_reflection_s5.png" width="45%"/></img>
<img src="OutputFiles\balls_shading_shadow_reflection_s5.png" width="45%"/></img>

___

### Acceleration  Data Structure

The Bounding Volume Hierarchy reduce the complexity of intersection check 
from linear to logarithmic. Due to the fact that 
intersection check is called recursively in the program, the acceleration data structure
greatly reduces the rendering time.


#### The rendering time of following images:

Shading + shadow + reflection

Teapot:     **433 seconds** (before) to **19 seconds** (after)

Balls:      **975 seconds** (before) to **22 seconds** (after)

<img src="OutputFiles\teapot_shading_shadow_reflection.png" width="45%"/></img>
<img src="OutputFiles\balls_shading_shadow_reflection.png" width="45%"/></img>

---

## Problems encountered:
There was a problem that bothered me for a long time when I was implementing 
the acceleration data structure. My function that constructs a BVH tree did not stop recursion.

It turns out that it is a problem with the Tracer constructor, 
which is provided in the framework.

Without the commented block of code I added later, the tracer will read the last surface
in the input file twice, resulting in two identical surfaces in the scene, 
and I cannot separate them into different child nodes because they have exactly
the same location.


```c++
Tracer::Tracer(const std::string &fname) {
    std::ifstream in(fname.c_str(), std::ios_base::in);
    std::string line;
    char ch;
    Fill fill;
    bool coloredlights = false;
    while (in) {
        /***************************/
        // if (in.eof()) break;
        /***************************/
        getline(in, line);
    }    
```

___

There was a picture rendered not realistically because of an error in 
calculating reflection, but it was a good-looking one.

<img src="OutputFiles\_error_balls_reflection.png" width="100%"/></img>

---

## Summary

This course helped me learned a lot about ray tracing as well as 
other knowledge of computer graphics. 

I will keep working on adding other features in my spare time, such as refraction or 
filters that can smooth the artifacts when rendering low-poly models.


