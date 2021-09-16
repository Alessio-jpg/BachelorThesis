# BachelorThesis

Introduction
--------------

This is a Python implementation of the Wavelet Scalar Quantization algorithm (WSQ). 

It handles encoding/decoding from/to all major image standards such as PNG, JPEG, TIP, BMP, etc . . .


Requirements
------------

* [Python3](https://www.python.org/downloads/),
* [numpy](https://numpy.org/install/),
* [scikit-image](https://scikit-image.org/docs/stable/install.html),
* [bitstring](https://pypi.org/project/bitstring/).

Usage:
----------

Encoding:

    python3 encoder.py [-d] <input> <output>
Decoding:

    python3 encoder.py [-d] <input> <output>

The option [-d] will print detailed informations during the entire encoding/decoding process

Examples:

    python3 encoder.py "fingerprint.png" "compressed_fingerprint.wsq"
    python3 decoder.py -d "compressed_fingerprint.wsq" "decoded_fingeerprint.png"


Test Images
-------------
The NIST image dataset useful for testing the WSQ specification can be found [here](https://www.nist.gov/programs-projects/wsq-certification-procedure).


