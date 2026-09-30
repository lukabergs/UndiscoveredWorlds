Can we convert these gradients/color ramps to something accurate if interpolating colors?

Review (2026-09-03):
These are already one-pixel PNG strips: temp has 496 samples and rain has 359. They preserve exact 8-bit RGB samples and work with the existing strip importer. Forward interpolation can reproduce the ramp, but physical values require an explicit minimum and increment: the images contain no temperature/rainfall calibration. Inverse colour matching currently selects the nearest RGB sample and cannot disambiguate repeated colours (src/io/map_imports.cpp).
