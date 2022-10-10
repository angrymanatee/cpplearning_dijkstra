# Overview of Dijkstra's Algorithm Implementation

I attempted to make this a bit more "object oriented" and less efficient as a learning exercise.  The `Edge` class is probably unnecessary, since I could have kept an `unsorted_map` of output nodes and the weight associated with each `Vertex` and ensured they were consistent in the initializer.  However, I wanted to play around with shared pointers a bit and make it feel more object oriented, where a user could traverse a graph by bouncing between the `Vertex` and `Edge` pointers.
 
I attempted to keep it reasonably efficient by passing references to objects and whatnot in the functions.  I read somewhere that modern C++ compilers have return variable optimization, so I was hoping the code could take advantage of this.  I haven't done any profiling to see if I succeeded.
 
Overall, it was a fun little project where I played around with shared pointers and classes.  I spent a lot of time looking into the std libraries, since it’s better to use built in functions than roll new code.  I also learned that it’s possible to run a function on an uninitialized object, and it can segfault when trying to access `this`.  That was a fun bug.  I return `-1` as a “not found” signal to the function caller, and the signal wasn’t handled and resulted in indexing `vector[-1]`.
