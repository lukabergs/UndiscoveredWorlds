Make different kinds of appearance modes, such as continuous, discrete... and let user choose palette.

Review (2026-09-03):
Editable map gradients already offer a Discrete checkbox, value/colour anchors, and .uws presets (src/app/rendering/map_appearance.cpp). Continuous interpolation is used when Discrete is off. Benchmark diagnostic exports still use fixed palettes; exposing these through the same controls remains a separate change.
