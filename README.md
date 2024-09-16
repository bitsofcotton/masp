# masp
make and apply internal states from/to input stream.

# Description
We make internal states for a type of automorphism of input stream.
Also we apply them into another input stream.

The map we make is in the bitsofcotton/p1 described form.
However, the +p, -p command only depends on 3 back input images.

Making internal states is enough one in bit stream meaning as in p1.

# Usage
    ./masp(32)?(mp)? [+-][ap] in0.ppm ...
    # we assume + for making internal states, - for applying the states.
    # we assume a only for whole image, p for 3-image back prediction.

# Whole graphics prediction with p1
We can learn masp... +a ((whole image)) &gt; learn.txt then, pp3n each pixel projection causes internal space correct prediction in p1 meaning. (we should get last quad image and revert them by learned vector then average them.)

However, this is internal states valid prediction, isn't the next one image out of the states meaning.

