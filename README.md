<!-- [![DOI](https://zenodo.org/badge/28657065.svg)](https://zenodo.org/badge/latestdoi/28657065) -->
<div id="top"></div>
<br />
<div align="center">
  <a href="https://github.com/lagoproject/anna">
    <img src="docs/images/lago-logo.png" alt="Logo" width="140">
  </a>
  <h3 align="center">The ANNA Framework</h3>
  <p align="center">
    A framework designed to analyze the data measured by the single Water Cherenkov Detectors (WCD) of the <a href="https://lagoproject.net/">LAGO</a> detection network. It follows the LAGO hierarchycal structured for the measured data.
    <br />
    <!-- <a href="https://github.com/lagoproject/anna"><strong>Explore the docs (soon) »</strong></a>
    <br /> -->
    <br />
    <a href="https://github.com/lagoproject/anna/issues">Request Feature</a>
    ·
    <a href="https://github.com/lagoproject/anna/issues">Report Bug</a>
    ·
    <a href="#Contact">Contact us</a>
</p>
</div>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
<br />
  <ol>
    <li><a href="#about-anna">About ANNA</a></li>
    <li><a href="#getting-started">Getting Started</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#proposed-features">Proposed Features</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## About ANNA
<!--
<style>
html, body {height: 100%;}
img {
height: auto;
width: auto;
}
img.rel {
  width: 33%;
}
</style>
-->

<!-- <img class="rel" src="./docs/images/flux-chacaltaya.png" alt="The seconday particle flux at Chacaltaya" width="200"><img class="rel" src="./docs/images/geomagnetic.png" alt="Charged particles trajectories" width="200"><img class="rel" src="./docs/images/wcd-muon.png" alt="A muon traversing the WCD " width="200">
-->

ANNA is a complete framework designed to analyse the data the signals produced by the secondary particles emerging from the interaction of singles, multiple and even the complete flux of primary cosmic rays with the atmosphere. These signals are measured for any particle detector located at a LAGO site.

ANNA is structured using different applications files and some headers where the main data classes are defined.

During its normal operation and the different maintenance modes of the <a href="https://lagoproject.net/">LAGO detection network</a>, LAGO produces datasets containing different types of data. All the LAGO datasets follow standardized schemes based on linked data as described in our <a href="https://lagoproject.github.io/DMP/">DMP</a>. 

<div id="measured"></div>

### Measured data (Ln)

The measured data LAGO dataset corresponds to any type of data that was measured in a LAGO site, including measurements of any type of radiation, atmospheric conditions, geomagnetic, solar irradiation, telemetry, and any other type of measurement produced by the LAGO detector and its peripherals. The LAGO measured data includes also any data product derived by any means of the previously described data. All the LAGO measured datasets are hierarchically tagged with the ***L<sub>n</sub>*** label, being ***n*** a sequential number, starting by *0*, that indicates the level of data processing: L<sub>0</sub> corresponds to the raw data, L<sub>1</sub> is the first level of analysed data, etc.

### ANNA main reference and citation

When using ANNA, please cite using the DOI included at the begininig of this file. 

H. Asorey for the LAGO Collaboration, _"The LAGO ANNA Data Analysis framework"_, [doi:_SOON_](https://doi.org/).

<p align="right">(<a href="#top">back to top</a>)</p>

## LAGO ANNA versions 

Currently, two electronic boards are consistenly used across the LAGO detection network. Two different versions were developed for [LAGO ACQUA](https://github.com/lagoproject/acqua), the LAGO data acquisition package. LAGO ANNA follows the same schema:

* LAGO ANNA v1.x.x: Compatible with LAGO ACQUA v1, for the electronic board based on the NEXYS-II FPGA. Current version is the last developed version, [LAGO ACQUA v1.5.0](https://github.com/lagoproject/acqua/releases/tag/acqua-v1.5.0). In general, LAGO ANNA v1 are compatible with data acquired using LAGO ACQUA v1. 
* LAGO ANNA v2.x.x: Compatible with LAGO ACQUA v2, designed for the electronic board based on the StemLab 125/14 RedPitaya electronic board. This is currently under actively development. 
* LAGO ANNA v3.x.x: This will be developed for the new LAGO EDGE detection schema. 

Major versions for both the LAGO ACQUA and the LAGO ANNA frameworks are developed in different branchs and identified with the corresponding tags. Master branch always points to the latest stable development. 


<!-- GETTING STARTED -->
## Getting Started

To get a local copy of the ANNA framework up and running follow these simple example steps.

### Prerequisites

#### System requirements

ANNA runs in any Linux based system, including those supported at iOS, raspberry-pis' and RedPitayas' linux-based OS. For Windows user, we strongly recommend to install some of your preferred linux distribution using virtualbox. 

The command for installing required packages depends on the OS architecture. 

In Fedora, Scientific Linux and CentOS, use 

`sudo yum install <package>`. 

In Ubuntu/Debian, including RPi and RPy SBC, use 

`sudo apt install <package>`. 

In both cases, you can also use the graphic package manager included in your preferred distro. 

ANNA requires the installation of a few and largerly common standard packages: 

* bash
* gcc
* make
* screen
* rsync
* git 

As a one-liner for Ubuntu/Debian:

```bash 
sudo apt install build-essential screen rsync git
```
(git is optional).

#### Dependencies

ANNA does not have any dependencies.

### Installation

1. If you are using git, just clone this repository:
   ```bash
   cd /path/to/ANNA/installation
   git clone https://github.com/lagoproject/anna.git
   ```
   Otherwise, you can also directly download ANNA without using git (in this case, you should reinstall ANNA for every upgrade):
   ```bash
   cd /path/to/ANNA/installation
   wget -c https://github.com/lagoproject/anna/archive/refs/heads/master.zip
   unzip master.zip
   rm master.zip
   ```
2. ANNA compiling is very simple:
   ```bash
   cd /path/to/ANNA/anna
   make
   ```
   1. During the first installation of ANNA (or if you need to install ANNA in a different directory), `make` will define the `$LAGO_ANNA` environment variable, that points to the ANNA current directory installation. Then, ANNA installer will add the definition of this variable to the user's local `.bashrc` and to the local `$PATH` environment variable:
   
      ```bash
      #
      ## Changes added by the ANNA suite on <installation date>
      #
      export ANNA="/path/to/ANNA/installation/ANNA"
      export LAGO_ANNA_VERSION="<LAGO ANNA current version>"
      export PATH="${ANNA}:$PATH"
      ```

If you follow the above described steps and everything works well, you should find some new executable files at the root `${LAGO_ANNA}` directory.


<p align="right">(<a href="#top">back to top</a>)</p>

### ANNA updates, releases, branchs and tags

ANNA is continously used, revised and updated within the [LAGO Collaboration](https://lagoproject.net). 

Unless you are a developer, we recommend to use only the latest ANNA release contained in the `master` branch of this repository. Stable versions are tagged and can be found in the [corresponding section of this repository](https://github.com/lagoproject/anna/tags).

Clone and install ANNA from `dev` or `dev-*` branches is strongly discouraged, as these branches are used for testing, bug correction and for the development of new features.

If you are using `git`, you can update ANNA just by doing: 

```bash
cd /path/to/ANNA/anna
git pull
make
```

Otherwise, you could just reinstall ANNA by following the [installation guide](#installation).

<p align="right">(<a href="#top">back to top</a>)</p>

## Usage

The ANNA framework follows the basic hierarchical structure of the <a href="#measured">LAGO Measured data</a>. It includes two basic applications, [`dump`](src/dump.cc) and [`example`](src/example.cc) intended as basic examples of the usage of the LAGO ANNA classes. 

All the ANNA applications have their own integrated help, accesible through the `-?` modifier. 

### General help and documentation

The documentation is currently under preparation and will be released soon. A brief description of the action of each code and the available options and modifiers can be seen by calling them with the `-?` modifier.

<!-- _For more examples, please refer to the [ANNA Documentation](docs)._ -->

<p align="right">(<a href="#top">back to top</a>)</p>

## Proposed features

See the [open issues](https://github.com/lagoproject/anna/issues) for a full list of proposed features (and known issues).

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- LICENSE -->
## License

ANNA is distributed under the [BSD-3 License](https://opensource.org/licenses/BSD-3-Clause). See the [LICENCE](LICENSE) for more information.

<p align="right">(<a href="#top">back to top</a>)</p>

<!-- CONTACT -->
## Contact

The ANNA framework is developed by the LAGO Collaboration. If you need to contact us, please complete our [contact form](https://lagoproject.net/contact.html). 

ANNA principal contact: [Dr Hernán Asorey (@asoreyh)](https://github.com/asoreyh)

Project Link: [https://github.com/lagoproject/anna](https://github.com/lagoproject/anna)

<p align="right">(<a href="#top">back to top</a>)</p>