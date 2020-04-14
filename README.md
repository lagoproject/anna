# ![The LAGO Project](http://lagoproject.net/images/lago-logo-90.png "The LAGO Project") The LAGO ANNA suite

| CODENAME			| ANNA  |
|-------------------|:------|
| COPYRIGHT			| (C) 2012-Today, The LAGO Project, [lagoproject.net](http://lagoproject.net)|
| LICENSE			| BSD-3-Clause |
| REPOSITORY		| https://github.com/lagoproject |
| CONTACT			| [lago@lagoproject.org](mailto:lago@lagoproject.org)|
| DESCRIPTION		| The LAGO project operates its WCD by using different types of electronics, but using a unified data structure. Here you will find the latest stable (prod) version of the data analysis suite of the LAGO Project |
| CONTRIBUTORS		| If you want to contribute to this project please [send us an email](mailto:lago@lagoproject.org)|


| CODE GUIDELINES	|		|
|-------------------|:------|
| FILE ENCODING		| UTF8 (please use <kbd>iconv -f YOUR_ENCODING -t UTF-8 file_to_convert > converted_file </kbd> before to push) |
| LANGUAGE			| English (preferred) |
| INDENT STYLE		| [Stroustrup](http://en.wikipedia.org/wiki/Indent_style#Variant:_Stroustrup) using 1 tab for 4 columns wide. check [here for vim setup](http://tedlogan.com/techblog3.html) |
|					| If you prefer, please use: <kbd>astyle -t4 -A4 -y file_to_convert</kbd> before to push
| VERSIONING		| Sequence-based identifiers, v<version>r<release>. First public release: **v1r0**
| INSTALL			| After installing dependences (see INSTALL), just *make*
| USAGE				| Please visit our [wikipage](http://wiki.lagoproject.org) (internal use only)|

The [Latin American Giant Observatory (LAGO)](http://lagoproject.org) is an extended Astroparticle Observatory at global scale. It is mainly oriented to basic research on three branches of Astroparticle physics: the Extreme Universe, Space Weather phenomena, and Atmospheric Radiation at ground level.


# Before to used it

## Libraries needed:

- sudo apt install openssh-server vim libreadline7 libreadline7-dev ntp xorg xorg-dev libusb-dev cc1111

- sudo apt install m4 automake autoconf gnu-standards autotools-dev autoconf-archive automake cvs subversion git mercurial manpages-dev gfortran build-essential gsl-bin libgd2-xpm-dev fort77 cmake rar libx11-dev rsync screen htop xclip dos2unix vpnc nmap pdfshuffler pdfsam

- sudo apt install python-numpy python-scipy python-matplotlib python-simpy python-simpy-gui imagemagick gnuplot-x11 astyle pdl

- sudo apt install hdf4-tools hdf5-tools libcfitsio-bin libcfitsio-dev libhdf5-dev libhdf4-dev pandoc
- sudo apt install python3-numpy python3-scipy ipython3 python3-matplotlib python3-simpy

