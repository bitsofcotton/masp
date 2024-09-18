# masp
make and apply internal states from/to input stream.

# Description
We make internal states for a type of automorphism of input stream.
Also we apply/invert them into another input stream.

Making internal states is enough one in bit stream meaning as in p1.

# Usage
    ./masp(32)?(mp)? [+-i] in0.ppm ...
    # we assume + for making internal states,
    #           - for applying into the internal states invariant.
    #           i for inverting image.

# Tips
    ./masp(32)?(mp)? + in0.ppm ... > L.txt
    ./masp(32)?(mp)? - another0.ppm < L.txt
    ...
    predg(32)?(mp)? another0.ppm-i4.ppm ...
    ./masp(32)?(mp)? i predg.ppm < L.txt
    #
    # This chain also
    #
    ./masp(32)?(mp)? + in0.ppm ... > L.txt
    python3 ../goki_check_cc/test.py ./masp(32)?(mp)? mprep L.txt another0.ppm ...
    predg(32)?(mp)? another0.ppm-i4.ppm ...
    ./masp(32)?(mp)? i predg.ppm < L.txt
    python3 ../goki_check_cc/test.py ./masp(32)?(mp)? mapply L.txt predg.ppm-i.ppm

This chain causes internal states depend predictions with goki_check_cc:test.py bit operation.

However, this is internal states valid prediction, isn't the next one image out of the states meaning.

Also, we can slide L.txt transision and predict them, however, we cannot apply this concern with another images.

