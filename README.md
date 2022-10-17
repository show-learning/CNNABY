# CNNABY
Build a CNN model with ABY framwork

We assume MPC(Secure Multi-Party Computation) is a solution to prove data security when using Cloud AI.
And we use ABY framwork to do that.

ABY framwork provides a security way to transform information.
If we make a CNN model by ABY, we can ensure the information we give to the model is safety. Don't need to worry about revealing.
It also mean we can build a CNN model AI by ABY, training or counting in a safety way.

If the AI model open, user can safety use the model, and the model builder can keep it training data.
Make people use the Al model without worried.

We choose the FSRCNN in YUV420 as target.

Using the layers and rewited code based on https://github.com/miladabd/DeepC

To prove ABY can be used in AI model.  

Demo version: Run about 16 days/1frame

Speedup version: Run about 1minuite/1frame

Fully version: Run about 1year/1frame ,but since the hardware limit(memory shortage) & time，we don't run the about

