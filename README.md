
# The SBL sandbox

The **sbl-sandbox** sources are a nice padded area in which we are safe
to play around with methods in _sparse Bayesian learning_. They are
essentially a set of lightweight solvers paired with some python
code that generates problem instances of the "_k_ spikes measured
with or without noise by an i.i.d. Gaussian random matrix" variety.

These sources support the following manuscript:

> Worley, B., _Scalable mean-field sparse Bayesian learning_,
> IEEE Transactions on Signal Processing, 2019, **in prep**.

## Usage

Bear with me, compilation is a bit non-standard.

To compile a source file, a problem instance, _e.g._ __src/instdef.hh__, is
paired with the core instance header __src/inst.hh__ and a solver,
_e.g._ __src/gs.cc__. Compiling would go a bit something like this:

```bash
g++ -std=c++14 -O3 -I. -I/path/to/eigen3 \
  -include src/instdef.hh -include src/inst.hh \
  src/gs.cc -o gs
```

Normally, manual compilation isn't necessary. The compilations needed
to perform all experiments are managed by __make.py__ using
[doit](http://pydoit.org).

## Preparation on *m5.metal*

The following commands are needed before running __make.py__ on EC2:
```bash
sudo yum -y groupinstall "Development Tools"
sudo amazon-linux-extras install python3
sudo pip3 install --upgrade pip
sudo pip3 install doit
```

Then, clone this repository and run:
```bash
git clone git@github.com:geekysuavo/sbl-sandbox.git
cd sbl-sandbox
./make.py -n 96
```

## Licensing

All sources in this repository are released under the
[MIT license](https://opensource.org/licenses/MIT). See the
[LICENSE.md](LICENSE.md) file for the complete license terms.

~ Brad.

