# Third-party notices

This file documents third-party material bundled with MitoSnake v0.9.0.

## Tiara models in `resources/tiara_models/`

MitoSnake ships classifier model files for use with
[Tiara](https://github.com/ibe-uw/tiara) (Karlicki, Antonowicz, and
Karnkowska, *Bioinformatics* 2022;38(2):344–350). `install.sh` copies them
into the `mitosnake_tiara` conda environment (unless
`--no-replace-tiara` is used). The copy step is implemented by
`scripts/replace_tiara_models.py`, which refers to the files as
"optimized" replacements for the models that ship with the Tiara package.

Layout (relative to `resources/tiara_models/`):

- `nnet/first_nnet_kmer_{4,5,6}.pkl`
- `nnet/second_nnet_kmer_{4,5,6,7}.pkl`
- `tfidf/k{4,5,6}-first-stage/{model.npy,params.txt}`
- `tfidf/k{4,5,6,7}-second-stage/{model.npy,params.txt}`

### Provenance (checked against ibe-uw/tiara `main`)

SHA-256 comparison with
https://github.com/ibe-uw/tiara (tree `main`, models under
`tiara/models/`):

- **First-stage TF-IDF** (`tfidf/k{4,5,6}-first-stage/model.npy` and
  `params.txt`) is **byte-identical** to the upstream defaults.
- **First-stage neural network** `nnet/first_nnet_kmer_4.pkl` is
  **byte-identical** to upstream. `first_nnet_kmer_5.pkl` and
  `first_nnet_kmer_6.pkl` have the same byte sizes as the upstream files
  (25,224,293 and 41,976,933 bytes); they were not re-hashed in this
  check but match the upstream names and sizes.
- **Second-stage TF-IDF and neural-network files are not identical** to
  that GitHub `main` snapshot. Example: upstream
  `k7-second-stage/params.txt` has `N:92336` and
  `data_names:data/mitochondria,data/plastids`; the bundled copy has
  `N:55044` and empty `data_names`. Second-stage `.pkl` files also differ
  in size from upstream.

Conclusion: first-stage files are the upstream Tiara defaults (confirmed
for TF-IDF and for `first_nnet_kmer_4.pkl`). Second-stage files are
**not** the current upstream defaults; they are the replacement set that
`install.sh` installs. This tree does not record the training command or
dataset used to produce the second-stage files, so they are described here
as **derived from Tiara** (same file names and layout; hashes differ).
They are **not** claimed to be a documented in-house retrain with a named
dataset.

All of these files remain subject to the Tiara MIT licence below.

## Tiara MIT License

The following is the licence text from
https://raw.githubusercontent.com/ibe-uw/tiara/master/LICENSE
(also present as `LICENSE` on branch `main`).

MIT License

Copyright (c) 2020 Institute of Evolutionary Biology, University of Warsaw

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
