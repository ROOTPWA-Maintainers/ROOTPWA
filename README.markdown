[TOC]


# ROOTPWA #

ROOTPWA is a toolkit for partial-wave analysis of multi-particle final
states produced in high-energy particle reactions. Based on
experimental and Monte-Carlo data that are read into the program in
the form of three-vector tuples and the set of partial waves specified
by the user in form of ASCII files, ROOTPWA determines the production
amplitudes of the various partial waves by performing an unbinned
extended maximum-likelihood fit that uses the full kinematic
information of the measured final-state particles and that takes into
account detector acceptance as well as reconstruction efficiency. For
the calculation of the partial-wave decay amplitudes the helicity
formalism is used. The code structure, however, is open for the
implementation of other formalisms.

Please note that ROOTPWA is still in development and currently in a
transition phase from an old deprecated code basis using a mixture of
C++, bash, and ROOT scripts to a new system that is based on C++ and
Python. The new system is intended to be faster, more versatile, more
flexible, and easier to maintain. The transition will also encompass
the migration of the main I/O data formats towards efficient and
self-describing ROOT data structures. We will try to minimize
incompatibilities and disruptions of workflows that will inevitably
occur along the way.

ROOTPWA is distributed under the [GPLv3](http://www.gnu.org/licenses/gpl.html).
For details please see the LICENSE file contained in this
package. ROOTPWA is hosted at SourceForge and can be found at

<http://sourceforge.net/projects/rootpwa/>

Below you will find instructions on how to install the external
libraries needed by ROOTPWA, how to get the ROOTPWA sources, how to
install the program, and how to contribute to the project.


***


# Installation Instructions #


## Prerequisites ##

ROOTPWA will usually run without problems on most (not too old) Linux
and MacOS X systems. Regular tests are performed on the latest Ubuntu
distribution and MacOS version as well as on the Scientfic Linux CERN
(SLC) distribution installed on the [_lxplus_](http://plus.web.cern.ch/plus/)
interactive logon service provided by CERN. Before you can build
ROOTPWA you need to install the following software packages.


### CMake ###

CMake is a cross-platform, open-source build system available from
<http://cmake.org>. In addition to the [reference documentation](http://www.cmake.org/cmake/help/documentation.html)
available at the project web page further information about CMake can
be found at

*   <http://www.cmake.org/cmake/help/cmake_tutorial.html>
*   <http://www.cmake.org/cmake/help/runningcmake.html>
*   <http://www.cmake.org/cmake/help/syntax.html>
*   <http://www.cmake.org/Wiki/CMake>
*   <http://rachid.koucha.free.fr/tech_corner/cmake_manual.html>
*   <http://mash-project.eu/wiki/index.php/CMake%3a_Getting_Started>

The minimum required CMake version is 2.8.0. In case your system
offers only outdated packages, you can quite easily compile CMake
yourself...

1.  Download the latest CMake release from
    <http://www.cmake.org/cmake/resources/software.html>

2.  Prepare the installation

    `> bootstrap --prefix=<install directory>`

    On multi-core machines you may speed up things by using the
    `--parallel=<# of processes>` option. If you plan to use the CMake
    GUI, you should set the `--qt-gui` flag. This of course requires
    that you have the QT libraries installed.

3.  Compile CMake

    `> make && make install`

4.  Add `<install directory>/bin` to your path or create an alias for cmake.

...or use the binary distribution package...

1.  Download the latest CMake binary (`.tar.gz`) from
    <http://www.cmake.org/cmake/resources/software.html> and unpack it to
    your preferred binary directory (e.g. `~/bin/`).

2.  Make sure the binary directory is in your path or create an alias for cmake.


### Boost ###

Part of the code relies on the Boost C++ template library which is
available at <http://www.boost.org>. Version 1.50.00 or higher is
required.  Boost is to a large extend a header-only library so that
usually nothing has to be build. Just install the respective Boost
packages for your platform or, in case you do not have administrator
privileges, extract the source archive to a location of your choice.

A more convenient way than downloading and extracting the tarball of
a certain Boost version is to clone the Boost git repository by running

    > git clone https://github.com/ryppl/boost-svn.git

and then to checkout the wanted release version via

    > cd boost-svn
    > git checkout Boost_1_50_0

This way one can easily switch to any existing Boost version and also
updates are much more convenient:

    > git checkout master
    > git pull
    > git checkout Boost_1_51_0

To find out the current version (a.k.a. branch tag) run

    > git describe --tags

The list of available versions (tags) is printed by

    > git tag


### ROOT ###

ROOT is an open-source data-analysis framework for high-energy and
nuclear physics and is available from <http://root.cern.ch>. Version
5.26 or higher is required and it must have been built with the
configure options

    --enable-mathmore
    --enable-minuit2

If `root-config --features` does _not_ list `mathmore` and `minuit2` you
need to re-configure and re-compile ROOT with these options.


### libconfig ###

We use the _libconfig_ config file parser written by Mark A. Lindner
available from <http://www.hyperrealm.com/libconfig/>. To our
knowledge there are unfortunately no pre-compiled packages for common
Linux platforms available, but compilation and installation are
straightforward. The easiest way is to install the library into the
same directory as the source

    > ./configure --prefix=/your/folder/to/install/libconfig
    > make && make install


### Python (optional) ###

In order to make scripting more powerful and flexible, it is planned
to Python-ify some of the ROOTPWA classes, so that they can be
interfaced directly. In the long term much of the house-keeping and
user-interface code that is currently scattered across several C++
programs, shell scripts, and ROOT scripts will be reimplemented in
Python.

The build system tries to find your Python installation
automatically. For this to work you need to have the `python`
executable in your path. In addition you also need to compile the
`Boost.Python` library (e.g. by running the supplied
`compileBoostLibraries.sh` script).  If the build system has found
your Python installation and the `Boost.Python` library, the Python
features are automatically enabled.

ROOTPWA requires Python 2.7. In case you do not have the possibility
to install the packages for your operating system, you may install
Python from source as outlined below:

1.  Download the source tarball from <http://www.python.org> and extract it to
    a directory of your choice.

2.  Compile Python.

    `> ./configure --enable-shared && make && make install && make test`

    Depending on whether you have administrator rights or not you
    might want to set the prefix accordingly (e.g. `` --prefix=`pwd` ``).
    In this case you also have to make sure to adapt the `PATH`
    and `LD_LIBRARY_PATH` environment variables accordingly.


### CUDA (optional) ###

If you have access to a CUDA capable nvidia Graphics Card (Shader Model 1.3 or
higher) you may want to install the CUDA framework (version 3.0 or
higher). The build system tries to find your CUDA installation based
on the `NVSDKCUDA_ROOT` environment variable. This should point to
`<your NVIDIA_GPU_Computing_SDK directory>/C`. If `NVSDKCUDA_ROOT` is
not set the build system assumes that the SDK is located in
`${HOME}/NVIDIA_GPU_Computing_SDK/C`. If a CUDA installation is found,
the CUDA features are enabled automatically.


### MPI (optional, experimental) ###

In order take advantage of the parallel nature of the computing
problems in PWA, it is planned to make some of the executables
MPI-aware, so that they run on multi-core machines as well as on MPI
PC-clusters. The build system tries to find your MPI installation
(openMPI recommended) automatically. In addition you also need to
compile the `Boost.MPI` libraries (e.g. by running the supplied
`compileBoostLibraries.sh` script). If the build system has found both
the MPI installation and the `Boost.MPI` libraries, the MPI features
are automatically enabled.


***


## Getting ROOTPWA ##

The ROOTPWA source code is available through the central [git
repository](https://sourceforge.net/p/rootpwa/code/) hosted at
SourceForge.  In order to get the sources you have to "clone" the
repository by running

    > git clone git://git.code.sf.net/p/rootpwa/code rootpwa-code

The command will download the code into the `rootpwa-code`
directory. This working copy is a git repository of its own which
contains _all_ past revisions of the project and which you may use for
code experiments. If you plan to contribute code, please read the
section "Contributing to ROOTPWA" below.

ROOTPWA is developed in multiple code branches, most of which contain
work in progress that might not work as expected. However, there are
dedicated stable branches that are intended to be used for real
analyses. Tested versions of these branches are identified by tags.

By default you will be in the "master" branch after you cloned the
repository. This branch is used for development so users are strongly
recommended to switch to a more stable branch. At the moment there is
one stable branch called `_v1`. You can get a list of tags for this
branch by running

    > git tag | grep <branch name>

The naming scheme for the tags is `<branch name>.<version number>`. In
order to get the revision that belongs to a certain tag run

    > git checkout <tag name>

If you are more adventurous and want to get the latest version of a
branch, run

    > git checkout -t origin/<branch name>


***


## Building ROOTPWA ##

Finally you are ready to build ROOTPWA from the sources.


### 1. Setup Build Environment ###

1.  Make sure that you have installed all the packages mentioned above
    for the correct architecture: If you compile on a 32 bit machine,
    _all_ external libraries need also to be compiled for 32
    bit. Similarly for 64 bit machines _all_ libraries (including
    ROOT!) have to be compiled for 64 bit. Mixed 32/64 bit compilation
    is not supported and will result in linker errors.

2.  Make sure that your ROOT environment is setup correctly. The build
    system utilizes the `root-config` executable to determine the ROOT
    environment and expects it to be in the path.

3.  Define the `LIBCONFIG` environment variable and point it to the
    libconfig installation directory. The build system expects the include
    files in `${LIBCONFIG}/include` and the libraries in
    `${LIBCONFIG}/lib`. `${LIBCONFIG}/lib` should be added to the
    `LD_LIBRARY_PATH` environment variable.

4.  If the Boost library was _not_ installed in the system directories
    you have to tell the build system, where it can find the Boost
    files by defining the `BOOST_ROOT` environment variable such that it
    points to the Boost top level directory (e.g. `export
    BOOST_ROOT=some/path/boost-svn`).

5.  Set the `ROOTPWA` environment variable to the path of the ROOTPWA top
    level directory. This usually is your git working copy.

ROOTPWA comes with some example `setup*.sh` scripts that define set
all the above mentioned environment variables and can be used as
templates. You should copy one of them and modify it according to your
environment.


### 2. Compiling ROOTPWA ###

1.  We use an out-of-source build strategy (one of the nice features of
    CMake). This means that all files created by the build system will
    reside in a separate directory outside of the directories that
    contain the source code. The build system expects this directory
    to be `${ROOTPWA}/build` which should already exist in your git
    working copy.

2.  The contents of the build directory _can be safely deleted_ whenever
    you feel like it, since it can be _always_ regenerated by (re)starting
    the build process. This also means that is not a good idea to put
    any valuable files into the `build` directory or one of its
    subdirectories. To start the build process run

    `> cd ${ROOTPWA}/build`

    `> cmake ..`
    
    `> make`

    Note that there is no `make install`.

    If you run on a multi-core machine, you may want to use `make -j`
    in order to run multiple compilation processes in parallel which
    usually speeds up the whole process significantly.

    You may influence the build process via four build options:

    1.  `RELEASE`: This is the default.
    2.  `DEBUG`: This generates debug symbols and switches off all optimizations.
    3.  `RELWITHDEBINFO`: This generates debug symbols but uses the same optimization level like `RELEASE`.
    4.  `MINSIZEREL`: This optimizes for minimum binary size.

    You can choose a different build option by using CMake's `-D` option, e.g.

    `> cmake -D CMAKE_BUILD_TYPE=DEBUG ..`

    `> cmake -D CMAKE_BUILD_TYPE=RELWITHDEBINFO ..`

    `> cmake -D CMAKE_BUILD_TYPE=MINSIZEREL ..`

    CMake will print the compiler flags used in the build process.
    For debugging purposes you can enable the CMake debug output by
    running

    `> make VERBOSE=1`

3.  In order to be able to run the ROOTPWA programs you need to add
    `${ROOTPWA}/build/lib` to your `LD_LIBRARY_PATH` environment variable. It is also a good idea to
    add `${ROOTPWA}/build/bin` to the path.

4.  If you have doxygen installed on your system you can build a
    html documentation by running

    `> make doxygen`

    The doxygen documentation will be generated in the directory
    `${ROOTPWA}/html-doc/html/`. It can be viewed with any web browser, e.g.

    `> firefox ${ROOTPWA}/html-doc/html/index.html`

    This doxygen documentation is also available at the SourceForge project site at
    <http://rootpwa.sourceforge.net/>. However, it is updated less frequently and might thus be slightly outdated.

    Registered SourceForge project members with the respective permissions can upload the
    documentation by doing

    `> cd ${ROOTPWA}/html-doc`

    `> sftp <username>,rootpwa@web.sourceforge.net`

    `sftp> cd htdocs`
   
    `sftp> put html/*`

    `sftp> exit`


***


# Contributing to ROOTPWA #

Contributions to the development of ROOTPWA are very welcome and can
be made in several ways:


## 1. Bug Reports / Feature Requests ##

If you find a problem with ROOTPWA or feel that a feature that you
would like to have is missing, you can file a [bug
report](https://sourceforge.net/p/rootpwa/bugs/) or a [feature
request](https://sourceforge.net/p/rootpwa/feature-requests/) on our
website. Before submitting a report, please check if your item has not
already been submitted by somebody else. Note that you can vote on
both bug reports and feature requests if you feel that their
resolution is of help to you.

If you already have a solution for a particular bug or have already
implemented a new feature, feel free to attach a patch to the
corresponding request. The ROOTPWA developers will review your changes
and consider them for implementation to the repository.


## 2. Contributing Code ##

If you would like to contribute to ROOTPWA to an extent which would be
impractical for patch files, you can fork the ROOTPWA repository. To
do this, you need to have a SourceForge account. While logged in, go
to the [ROOTPWA repository](https://sourceforge.net/p/rootpwa/code/)
and click "Fork" on the left. This will give you your own copy of the
ROOTPWA git repository on SourceForge which you can use to implement
your contribution. As soon as you judge your work ready for
integration, you can issue a merge request by clicking "Request Merge"
in your copy of the repository. The ROOTPWA developers will then
consider your changes for implementation.

In case you are new to git and need information on how to work with
it, we recommend the [git book](http://git-scm.com/book), which gives
a thorough and comprehensible introduction to git.

Some information which might be of importance:


#### 1. Branch Layout ####

The ROOTPWA repository uses the master branch for development, meaning
that the latest development version is generally in the master
branch. For older versions, there are maintenance branches name "_v1",
"_v2", etc.. Versions which are deemed stable are tagged with a
version number, e.g. <branch name>.<version number>
(e.g. _v1.12). Most work is to be done in topic branches whose name
should _not_ start with an `_`. If you are unsure which branch your
work should be based upon, feel free to contact the ROOTPWA developers
for support.


#### 2. Commit Message Format ####

ROOTPWA follows a policy for the commit messages, which is close to
the standard git commit message policy. It consists of the following
rules:

- The commit message must consist of at least three lines
- The first line of the commit message must not be longer than 80 characters.
- The second line of the commit message must be empty
- All remaining lines (of which there must be at least one) must not be longer than 80 characters.

Please note that the ROOTPWA git repository uses a server hook to
enforce this policy. Should any of the commits you would like to see
integrated not fullfil these requirements, all of them will be rejected
and you will have to re-edit the problematic messages. To reduce
problems, there is a git hook in the repository which you can use on
the client side to check all the commit messages already when
commiting. You can install it by running

    > cd $ROOTPWA/gitUtils
    > ./installClientHook.sh

Please note that the `ROOTPWA` environment variable has to be set for
the script to work.
