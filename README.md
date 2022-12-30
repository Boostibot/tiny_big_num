# Tiny Big Num
*Big num library paricularly optimized for small big numbers*

Tiny big num is a single header, standard only library with minimal STL dependencies. **Currently in developement**

## No memory allocations
Storage is provided to each function up front and the required amount can be queried via the appropriate constant time function. This is level of control particularly important for preventing cache evictions.

## Trivial API types
Makes no assumptions about the container used for storing numbers and required no additional meta data keeping. Numbers are represented by a pointer to data and size. These are usually kept together in a triviall POD struct called Slice.

## Cache friendly
Significant effort was put into reducing the amount of active memory usage especially copying of temporary results. Operations are preferably performed in the lowest possible amount of passes through the data and all basic algorithms can be performed in place. 
Examples of this include: 
- Quadratic multiplication (which is used as a last step of karatsuba multiplication) performs the multiplication and summation of the temporary in a single pass. This saves one pass through most of the output per digit of multiplied number (see fused_mul_quadratic). We use the same trick in few other functions. 
- Karatsuba multiplication uses the knowledge of output buffers to place two of the three summands directly into place requiring only one addition of the three and 2/3 less storage.
 - Exponentiation by squaring algorithm computes exactly how many iterations will be necessary. We than use this to determine into where to place the first result so that the final result will be in the output storage and no copies will be necessary.

## Extremely strict
Each algorithm starts by asserting exactly what is allowed including possible data aliasing. Corner cases are typically solved by not allowing them. While this results in harder to comform to api it removes the need for any internal validity checking and possible error code checking. It also means that if your code didnt crash due to an assert it will complete successfully.




