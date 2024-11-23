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
    # also we assume in0.ppm ... is goki_check_cc:test.py bit generated images.

# Tips
    ./masp(32)?(mp)? + in0.ppm ... > L.txt
    ./masp(32)?(mp)? - another0.ppm ... < L.txt
    ddpmopt(32)?(mp)? 0 another0.ppm-4.ppm ...
    ./masp(32)?(mp)? i <another0.ppm-height> predg.ppm ... < L.txt

This chain causes internal states depend predictions with goki_check_cc:test.py bit operation.

However, this is internal states valid prediction, isn't the next one image out of the states meaning.

Also, we can slide L.txt transision and predict them, however, we cannot apply this concern with another images.

# Tips (2)
If we make L.txt by partial one and grows L2.txt, ..., we can predict next one step by matrix predictions.

Either, when the size we use for the L.txt:matrix.cols() grows up, the prediction causes larger dimension specific values on the stream, we don't know what the result means.

