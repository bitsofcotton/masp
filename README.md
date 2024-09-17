# masp
make and apply internal states from/to input stream.

# Description
We make internal states for a type of automorphism of input stream.
Also we apply them into another input stream.

The map we make is in the bitsofcotton/p1 described form.
However, the +p, -p command only depends on 3 back input images.

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
    ddpmopt(32)?(mp)? another0.ppm-i4.ppm ...
    ./masp(32)?(mp)? i predg.ppm < L.txt

This chain causes internal states depend predictions.

However, this is internal states valid prediction, isn't the next one image out of the states meaning.

