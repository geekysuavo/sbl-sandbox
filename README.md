
# The SBL sandbox

Fix me...

These sources support the following manuscript, _in preparation for
submission to_:

> Worley, B., _Scalable mean-field sparse Bayesian learning_,
> IEEE Transactions on Signal Processing, 2019.

## Usage

Bear with me, compilation is a bit non-standard.

To compile a source file, a problem instance, _e.g._ __src/instdef.hh__, is
paired with the core instance header __src/inst.hh__ and a solver,
_e.g._ __src/gs.cc__. Compiling would go a bit something like this:

```bash
g++ -std=c++14 -I. -I/path/to/eigen3 \
  -include src/instdef.hh -include src/inst.hh \
  src/gs.cc -o gs
```

Normally, this compilation isn't necessary. The compilations needed
to perform all experiments are managed by __make.py__ using
[doit](http://pydoit.org).

## Licensing

All sources in this repository are released under the
[MIT license](https://opensource.org/licenses/MIT). See the
[LICENSE.md](LICENSE.md) file for the complete license terms.

~ Brad.

