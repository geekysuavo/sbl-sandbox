
# EP-SBL

Fix me...

These sources support the following manuscript, _in preparation for
submission to_:

> Worley, B., _Expectation Propagation for sparse Bayesian learning_,
> Journal of Machine Learning Research, 2018.

## Introduction

Fix me...

### Compilation

Bear with me, compilation is a bit non-standard.

To compile a source file, a problem instance, _e.g._ __instance.hh__, is
paired with a solver, _e.g._ __gs.cc__. To build a solver against the
default instance (__instance.hh__), you can run:

```bash
./build.py solver gs
```
which will create an executable file __gs__ in the current working
directory. To build a solver for another instance, you must specify
the instance identifier, for example:

```bash
./build.py solver ep AZ5437
```
which uses the instance file __instanceAZ5437.hh__ and creates an
executable file __epAZ5437__ in the current working directory.
Instance files may be created using __build.py__ like so:

```bash
./build.py instance k=13 > instanceK13.hh
```
which overrides the default _k_ parameter in __instance.hh__.

## Licensing

All sources in this repository are released under the
[MIT license](https://opensource.org/licenses/MIT). See the
[LICENSE.md](LICENSE.md) file for the complete license terms.

~ Brad.

